/*************************************************************************
 * Copyright (C) 2017,  Alexei Sibidanov                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or	 *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,	 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of	 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	 *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License	 *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *************************************************************************/

/*************************************************************************
 * The implementation of the Linear Congruential Random Number           *
 * Generator with a large modulus defined by the recurrence:             *
 * x_{i+1} = x_{i} * A mod m                                             *
 * where the recurrence parameters are based on the RANLUX generator:    *
 * the base b = 2^24                                                     *
 * the modulus m = b^24 - b^10 + 1 = 2^576 - 2^240 + 1                   *
 * the multiplier A is a power of a -- A = a^p mod m,                    *
 * where a = m - (m-1)/b                                                 *
 *         = b^24  - b^23  - b^10  + b^9   + 1                           *
 *         = 2^576 - 2^552 - 2^240 + 2^216 + 1                           *
 * is the multiplicative inverse of the base b = 2^24                    *
 * i.e. a * b mod m = 1                                                  *
 *************************************************************************/
#include <stdint.h>

#pragma once

#define   likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

class ranluxpp {
protected:
  uint64_t _x[9]; // state vector
  uint64_t _A[9]; // multiplier
  // The original choice when constructing random doubles is to use 
  // 52 of a possible 53 random bits per double,
  // which results in 11 random doubles per state vector.
  // For this option set:
  // _n_packed_doubles = 11
  // Alternatively one can use the full 53 random bits per double,
  // but this only gives 10 random doubles per state vector.
  // For this option set:
  // _n_packed_doubles = 10
  // NB other values of _n_packed_doubles are not valid and will not compile.
  static const int _n_packed_doubles = 11; //can be either 11 (default) or 10
  uint64_t _doubles[_n_packed_doubles]; // cache for double precision numbers 
  uint32_t _floats[24];  // cache for single precision numbers 
  uint32_t _dpos; // position in cache for doubles
  uint32_t _fpos; // position in cache for floats

  // get a = m - (m-1)/b = 2^576 - 2^552 - 2^240 + 2^216 + 1
  static const uint64_t *geta();

  // fill the cache with float type numbers
  void nextfloats();
  
  // fill the cache with double type numbers
  void nextdoubles();

  // transfrom the binary state vector of LCG to 24 floats 
  void unpackfloats(float *a);

  // transfoom the binary state vector of LCG
  // to n doubles (where n = 10 or 11)
  template <int n> void unpackdoubles(double *d);
  
public:
  // The LCG constructor:
  // seed -- jump to the state x_seed = x_0 * A^(2^96 * seed) mod m
  // p    -- the exponent of to get the multiplier A = a^p mod m 
  ranluxpp(uint64_t seed, uint64_t p);
  ranluxpp(uint64_t seed) : ranluxpp(seed, 2048){}

  // get access to the state vector
  uint64_t *getstate() { return _x;}

  // get access to the multiplier
  uint64_t *getmultiplier() { return _A;}

  // seed the generator by
  // jumping to the state x_seed = x_0 * A^(2^96 * seed) mod m
  // the scheme guarantees non-colliding sequences 
  void init(unsigned long int seed);

  // set the multiplier A to A = a^2048 + 13, a primitive element modulo
  // m = 2^576 - 2^240 + 1 to provide the full period (m-1) of the sequence.
  void primitive();

  // produce next state by the modular mulitplication
  void nextstate();

  // return single precision random numbers uniformly distributed in [0,1).
  float operator()(float __attribute__((unused))) __attribute__((noinline)){
    if(unlikely(_fpos >= 24)) nextfloats();
    return *(float*)(_floats + _fpos++);
  }
  
  // return double precision random numbers uniformly distributed in [0,1).
  double operator()(double __attribute__((unused))) __attribute__((noinline)){
    if(unlikely(_dpos >= _n_packed_doubles)) nextdoubles();
    return *(double*)(_doubles + _dpos++);
  }

  // Fill array size of n by single precision random numbers uniformly
  // distributed in [0,1).
  void getarray(int n, float *a);

  // Fill array size of n by double precision random numbers uniformly
  // distributed in [0,1).
  void getarray(int n, double *a);

};

// unpack state into 11 double precision numbers
// 52 bits out of possible 53 bits are random
template <> inline void ranluxpp::unpackdoubles<11>(double *d) {
  const uint64_t
  one = 0x3ff0000000000000, // exponent
  m   = 0x000fffffffffffff; // mantissa
  uint64_t *id = (uint64_t*)d;
  id[ 0] = one | (m & _x[0]);
  id[ 1] = one | (m & ((_x[0]>>52)|(_x[1]<<12)));
  id[ 2] = one | (m & ((_x[1]>>40)|(_x[2]<<24)));
  id[ 3] = one | (m & ((_x[2]>>28)|(_x[3]<<36)));
  id[ 4] = one | (m & ((_x[3]>>16)|(_x[4]<<48)));
  id[ 5] = one | (m & ( _x[4]>> 4));
  id[ 6] = one | (m & ((_x[4]>>56)|(_x[5]<< 8)));
  id[ 7] = one | (m & ((_x[5]>>44)|(_x[6]<<20)));
  id[ 8] = one | (m & ((_x[6]>>32)|(_x[7]<<32)));
  id[ 9] = one | (m & ((_x[7]>>20)|(_x[8]<<44)));
  id[10] = one | (m & ( _x[8]>> 8));

  for(int j=0;j<11;j++) d[j] -= 1;
}

// unpack state into 10 double precision numbers
// 53 bits out of possible 53 bits are random
template <> inline void ranluxpp::unpackdoubles<10>(double *d) {
  // mask to select 53 right-most bits, i.e. integer in range [0,2^53)
  const uint64_t m = 0x001fffffffffffff;
  //  2^-53 in double precision
  const double sc = 1. / (UINT64_C(1) << 53);
  // construct integer in range [0,2^53) using 53 random bits
  // then multiply by 2^-53 to get double in range [0,1)
  d[0] = sc * (m & ((_x[0]>>0)));
  d[1] = sc * (m & ((_x[0]>>53)|(_x[1]<<11)));
  d[2] = sc * (m & ((_x[1]>>42)|(_x[2]<<22)));
  d[3] = sc * (m & ((_x[2]>>31)|(_x[3]<<33)));
  d[4] = sc * (m & ((_x[3]>>20)|(_x[4]<<44)));
  d[5] = sc * (m & ((_x[4]>>9)));
  d[6] = sc * (m & ((_x[4]>>62)|(_x[5]<<2)));
  d[7] = sc * (m & ((_x[5]>>51)|(_x[6]<<13)));
  d[8] = sc * (m & ((_x[6]>>40)|(_x[7]<<24)));
  d[9] = sc * (m & ((_x[7]>>29)|(_x[8]<<35)));
}
