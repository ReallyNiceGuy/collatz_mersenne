#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <csignal>
#include <boost/io/ios_state.hpp>
#include <ios>
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
  bool is_odd() const { return num[0]&1;}
  bool is_one() const { return num.size()==1 && num[0]==1;}
};

ostream& write(ostream& ofs, const bignum& n)
{
  vector<uint64_t>::size_type len{n.num.size()};
  vector<uint64_t>::size_type cap{n.num.capacity()};
  ofs.write((char*)&cap,sizeof(cap));
  ofs.write((char*)&len,sizeof(len));
  ofs.write((char*)&n.num[0],sizeof(uint64_t)*len);
  return ofs;
}

istream& read(istream& ifs, bignum& n)
{
  vector<uint64_t>::size_type len{};
  vector<uint64_t>::size_type cap{};
  ifs.read((char*)&cap,sizeof(cap));
  ifs.read((char*)&len,sizeof(len));
  n.num.clear();
  n.num.reserve(cap);
  n.num.resize(len);
  ifs.read((char*)&n.num[0],sizeof(uint64_t)*len);
  return ifs;
}

void save(const std::string& fn, const bignum&n, double elapsed, uint64_t count, uint64_t zero_run)
{
  std::ofstream ofs(fn, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
  write(ofs,n);
  ofs.write((char*)&count,sizeof(count));
  ofs.write((char*)&elapsed,sizeof(elapsed));
  ofs.write((char*)&zero_run,sizeof(zero_run));
  std::cerr << "Elapsed: " << elapsed << ", steps: " << count << ", bits (approx): " << n.num.size()*64 << std::endl;
}

bool load(const std::string& fn, bignum&n, boost::chrono::system_clock::time_point&t, uint64_t& count, uint64_t& zero_run)
{
  std::ifstream ifs(fn, std::ios_base::binary | std::ios_base::in);
  if (ifs.fail()) return false;
  read(ifs,n);
  if (ifs.fail()) return false;
  ifs.read((char*)&count,sizeof(count));
  if (ifs.fail()) return false;
  double d;
  ifs.read((char*)&d,sizeof(d));
  if (ifs.fail()) return false;
  t = boost::chrono::system_clock::now() -
      boost::chrono::duration_cast<boost::chrono::system_clock::time_point::duration>(boost::chrono::duration<double>(d));
  ifs.read((char*)&zero_run,sizeof(zero_run));
  return !ifs.fail();
}

std::ostream& operator<<(std::ostream& o, const bignum& n)
{
  boost::io::ios_flags_saver  ifs( o );
  for(int i=n.num.size();i>0;--i)
  {
    o << " " << std::hex << std::setfill('0') << std::setw(16) << n.num[i-1];
  }
  return o;
}

uint64_t bignum::by2n(uint64_t& zero_run)
{
  uint64_t i{};
  if (num[0] == 0) //full item right shift, cheap
  {
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

  int lshift = __builtin_ctzll(num[0]); // How many zeroes left?
  int rshift = 64-lshift;
  auto last = num.size()-1;
  for(size_t i=0;i<last;++i)
  {
    num[i]=(num[i]>>lshift)|(num[i+1]<<rshift);
  }
  num[last]>>=lshift;
  if (num[last]==0 && num.size()>1) num.pop_back(); //remove MSI if 0
  auto steps=lshift+i;
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
  if (num[last]==0 && num.size()>1) num.pop_back(); //remove MSI if 0
  return 2;
}

uint64_t collatz(bignum& n, uint64_t steps, uint64_t& zero_run)
{
  while (!(n.is_one() || interrupted))
  {
    if (n.is_odd()) steps+=n.x3p1by2();
    else steps+=n.by2n(zero_run);
  }
  return steps;
}

bignum mersenne(int power)
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
    if (pos == argv[0])
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
      n=mersenne(val);
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
      std::cout << val << "," << steps <<"," << "\"" << min << "m" << sec << "s\"," << dur.count() << ", 2^" << zero_run << std::endl;
      std::remove(cache.c_str());
    }
  }
}
