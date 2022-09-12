#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <boost/io/ios_state.hpp>
#include <ios>

using std::vector;
using std::size_t;
typedef unsigned __int128 uint128_t;

struct bignum
{
  vector<uint64_t> num;
  // (3*x + 1) / 2
  int x3p1by2();
  // x / (2^n)
  int by2n();
  bool is_odd() const { return num[0]&1;}
  bool is_one() const { return num.size()==1 && num[0]==1;}
};

std::ostream& operator<<(std::ostream& o, const bignum& n)
{
  boost::io::ios_flags_saver  ifs( o );
  for(int i=n.num.size();i>0;--i)
  {
    o << " " << std::hex << std::setfill('0') << std::setw(16) << n.num[i-1];
  }
  return o;
}

int bignum::by2n()
{
  int i{};
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
  for(size_t i=0;i<num.size()-1;++i)
  {
    num[i]=(num[i]>>lshift)|(num[i+1]<<rshift);
  }
  num.back()>>=lshift;
  if (num.size()>1 && num[num.size()-1]==0 ) num.pop_back();
  return lshift+i;
}

int bignum::x3p1by2()
{
  //odd * 3 + 1 will always be even, so
  //avoid work and do 1 right shift at the same time
  //saves a couple of percent
  if (num.back()>0x4000000000000000ULL) num.push_back(0);
  uint128_t res{((uint128_t)3)*num[0]+1}; //3x+1 in uint128_t, we save the carry on the upper part
  num[0] = res; //lower part of 3x+carry
  res>>=64; //pull carry into lower uint64_t
  for(size_t i=1;i<num.size();++i)
  {
    res=((uint128_t)3)*num[i]+res; //3x+carry in uint128_t, we save the carry on the upper part
    num[i] = res; //lower part of 3x+carry
    res>>=64; //pull carry into lower uint64_t
    num[i-1]=(num[i-1]>>1)|(num[i]<<63); //divide previous by 2 and pull lowest bit from this one
  }
  num.back()>>=1; //divide MSI by 2
  if (num.size()>1 && num[num.size()-1]==0 ) num.pop_back(); //remove MSI if 0
  return 2;
}

uint64_t collatz(bignum& n)
{
  uint64_t steps{};
  while (!n.is_one())
  {
    if (n.is_odd()) steps+=n.x3p1by2();
    else steps+=n.by2n();
  }
  return steps;
}

bignum mersenne(int power)
{
  int rest=(power)%64;
  int items=(power)/64;
  bignum ret;
  if (items)
  {
    ret.num = vector<uint64_t>(items,0xFFFFFFFFFFFFFFFFULL);
  }
  if (rest)
  {
    ret.num.push_back(0);
    for(int i=0;i<rest;++i) ret.num.back()|=(1ULL<<i);
  }
  return ret;
}

int main(int argc, char **argv)
{
  if (argc > 1)
  {
    char *pos;
    int val = strtol(argv[1], &pos, 10);
    if (pos == argv[0])
    {
      std::cerr << "Needs a number as parameter" << std::endl;
      exit(1);
    }
    bignum n{ mersenne(val) };
    auto ret{collatz(n)};
    std::cout << ret << std::endl;
  }
}