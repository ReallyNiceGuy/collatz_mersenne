#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <csignal>
#include <ios>
#include <cctype>
#include <boost/chrono.hpp>

volatile int interrupted{0};

void signal_handler(int signal)
{
  interrupted = signal;
}

using std::vector;
using std::size_t;
using std::ostream;
using std::istream;
typedef unsigned __int128 uint128_t;

struct bignum
{
  vector<uint64_t> num;
  // (3*x + 1) / 2
  uint64_t x3p1by2();
  // x / (2^n)
  uint64_t by2n(uint64_t &zero_run);
  void rshift(int count);
  bool is_odd() const { return num[0]&1;}
  bool is_one() const { return num.size()==1 && num[0]==1;}
  bool is_zero() const { return num.size()==1 && num[0]==0;}
  static bignum mersenne(int power);
  static bignum two_np1(int power);
};

std::ostream& operator<<(std::ostream& o, const bignum& n);
std::istream& operator>>(std::istream& is, bignum& n);

void save(const std::string& fn, const bignum&n, double elapsed, uint64_t count, uint64_t zero_run)
{
  std::ofstream ofs(fn, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
  ofs << n << ' '
      << count << ' '
      << elapsed << ' '
      << zero_run;
  std::cerr << "Elapsed: " << elapsed << ", steps: " << count << ", bits (approx): " << n.num.size()*64 << std::endl;
}

bool load(const std::string& fn, bignum&n, boost::chrono::system_clock::time_point&t, uint64_t& count, uint64_t& zero_run)
{
  std::ifstream ifs(fn, std::ios_base::binary | std::ios_base::in);
  if (ifs.fail()) return false;
  ifs >> n;
  if (ifs.fail()) return false;
  ifs >> count;
  if (ifs.fail()) return false;
  double d;
  ifs >> d;
  if (ifs.fail()) return false;
  t = boost::chrono::system_clock::now() -
      boost::chrono::duration_cast<boost::chrono::system_clock::time_point::duration>(boost::chrono::duration<double>(d));
  ifs >> zero_run;
  return !ifs.fail();
}

std::ostream& operator<<(std::ostream& o, const bignum& n)
{
  std::ios_base::fmtflags f( o.flags() );
  o << std::hex;
  for(auto digit=n.num.rbegin(); digit!=n.num.rend();++digit)
  {
    o << std::setfill('0') << std::setw(16) << (*digit);
  }
  o.flags( f );
  return o;
}

std::istream& operator>>(std::istream& is, bignum& n)
{
  std::vector<char> buffer;
  std::istream::sentry s(is);
  if (s) while (is.good()) {
    int c = is.get();
    if (std::isdigit(c) ||
        (c>='a' && c<='f') ||
        (c>='A' && c<='F')) buffer.push_back(c);
    else
    {
      is.unget();
      break;
    }
  }
  if (buffer.size()==0)
  {
    is.setstate(std::ios::failbit);
  }
  else
  {
    n.num.clear();
    auto len = buffer.size();
    auto items = len/16;
    auto rest = len%16;
    buffer.push_back(0);
    for(decltype(items) i=0;i<items;++i)
    {
      auto pos=len-((i+1)*16);
      n.num.push_back(strtoull(&buffer[pos],nullptr,16));
      buffer[pos]=0;
    }
    if (rest)
    {
      n.num.push_back(strtoull(&buffer[0],nullptr,16));
    }
  }
  return is;
}

void bignum::rshift(int count)
{
  if (count==0) return;
  int lshift = 64-count;
  auto last = num.size()-1;
  for(size_t i=0;i<last;++i)
  {
    num[i]=(num[i]>>count)|(num[i+1]<<lshift);
  }
  num[last]>>=count;
  if (num[last]==0) num.pop_back(); //remove MSI if 0, LSI can never be 0
}

uint64_t bignum::by2n(uint64_t& zero_run)
{
  uint64_t i{};
  if (num[0] == 0) //full item right shift, cheap
  {
    if (num.size() == 1) return 0; //Nothing to do, number is 0
    //This is very unlikely to happen at all, but as
    // __builtin_ctzll is undefined for 0, I might as
    // well use the guard check of 0 for something useful

    //As we are here already, we can at least try to find as
    // many zero items as possible in sequence
    auto end = num.begin();
    do {
      ++end;
      i+=64;
    }
    while ( end != num.end() && *end == 0);
    num.erase(num.begin(), end); //and remove them in one full swoop
  }

  int count = __builtin_ctzll(num[0]); // How many zeroes left?
  rshift(count);
  auto steps=count+i;
  if (steps>=zero_run) zero_run = steps+1;
  return steps;
}

uint64_t bignum::x3p1by2()
{
  //odd * 3 + 1 will always be even, so
  //avoid work and do 1 right shift at the same time
  //saves a couple of percent
  if (num.back()>0x5555555555555554ULL) num.push_back(0);
  auto last = num.size()-1;
  uint128_t res{((uint128_t)3)*num[0]+1}; //3x+1 in uint128_t, we save the carry on the upper part
  num[0] = res; //lower part of 3x+carry
  res>>=64; //pull carry into lower uint64_t
  for(size_t i=1;i<num.size();++i)
  {
    res+=((uint128_t)3)*num[i]; //3x+carry in uint128_t, we save the carry on the upper part
    num[i] = res; //lower part of 3x+carry
    res>>=64; //pull carry into lower uint64_t
    num[i-1]=(num[i-1]>>1)|(num[i]<<63); //divide previous by 2 and pull lowest bit from this one
  }
  num[last]>>=1; //divide MSI by 2
  if (num[last]==0) num.pop_back(); //remove MSI if 0. LSI can never be 0
  return 2;
}

uint64_t collatz(bignum& n, uint64_t steps, uint64_t& zero_run)
{
  if (n.is_zero()) return 0;
  while (!(n.is_one() || interrupted))
  {
    steps += n.is_odd() ? n.x3p1by2() : n.by2n(zero_run);
  }
  return steps;
}

bignum bignum::mersenne(int power)
{
  if (power == 0)
    return bignum{.num{0}};
  int rest=(power)%64;
  int items=(power)/64;
  bignum ret;
  if (items)
  {
    ret.num = vector<uint64_t>(items,0xFFFFFFFFFFFFFFFFULL);
  }
  if (rest)
  {
    ret.num.push_back((1ULL<<rest)-1);
  }
  std::cerr << "bignum is aligned: " << !(((intptr_t)ret.num.data())%8) << std::endl;
  return ret;
}

bignum bignum::two_np1(int power)
{
  int rest=(power)%64;
  int items=(power)/64;
  bignum ret;
  ret.num = vector<uint64_t>(items+1);
  ret.num.back()=1ULL<<rest;
  ret.num[0]+=1;
  std::cerr << "bignum is aligned: " << !(((intptr_t)ret.num.data())%8) << std::endl;
  return ret;
}

int main(int argc, char **argv)
{
  if (argc > 1)
  {
      // Install a signal handler
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
    std::signal(SIGHUP, signal_handler);
    char *pos;
    int val = strtol(argv[1], &pos, 10);
    uint64_t steps{};
    uint64_t zero_run{};
    if (pos == argv[1])
    {
      std::cerr << "Needs a positive number as parameter" << std::endl;
      exit(1);
    }
    if (val<=0)
    {
      std::cerr << "Needs a positive number as parameter" << std::endl;
      exit(1);
    }
    std::string cache{argv[1]};
    cache+=".cache";
    boost::chrono::system_clock::time_point start;
    bignum n;
    if (!load(cache,n,start, steps, zero_run))
    {
      n=bignum::mersenne(val);
      std::cerr << "starting from scratch: " << argv[1] << std::endl;
      start = boost::chrono::system_clock::now();
    }
    else
    {
      std::cerr << "loaded cache file: " << cache << std::endl;
    }

calculate:
    steps = collatz(n, steps, zero_run);
    if (interrupted)
    {
      boost::chrono::duration<double> dur = boost::chrono::system_clock::now() - start;
      std::cerr << "\ninterrupted, saving cache file: " << cache << std::endl;
      save(cache,n, dur.count(), steps, zero_run);
      if (interrupted == SIGHUP)
      {
        interrupted = 0;
        goto calculate;
      }
    }
    else
    {
      boost::chrono::duration<double> dur = boost::chrono::system_clock::now() - start;
      auto sec = dur.count();
      int min = sec/60;
      sec = sec-(min*60);
      std::cout << val << "," << steps <<"," << "\"" << min << "m" << sec << "s\"," << dur.count() << "," << zero_run << std::endl;
      std::remove(cache.c_str());
    }
  }
}
