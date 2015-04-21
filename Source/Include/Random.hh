// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/Random.hh
//
// This header file contains the declarations for the MT19937 object and its
// wrapper class ParallelMT

#ifndef _RANDOM_HH_
#define _RANDOM_HH_

#include <vector>

#define NN       312
#define MM       156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM       0xFFFFFFFF80000000ULL //  most significant 33 bits
#define LM       0x7FFFFFFFULL         // least significant 31 bits

namespace Gaia {

// This is a 64-bit version of Mersenne Twister pseudorandom number
// generator. Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura.
// Please see the acknowledgement of use at the bottom of this file!
class MT19937 {

public:

	MT19937(){}
	// construct via single seed value
	MT19937(unsigned long long seed = 19650218ULL);
	// construct via array of seed values
	MT19937(unsigned long long init_key[], unsigned long long key_length);
	~MT19937(){};

	unsigned long long RandomInteger();
	double RandomReal();

protected:

	// the array for the state vector
	unsigned long long mt[NN];
	// == NN+1 means mt[NN] uninitialized
	int mti;
};

// the following class manages an array of MT19937 objects
class ParallelMT {

public:

	ParallelMT(){}
	ParallelMT(const int threads = 1,
        const unsigned long long first_seed = 19650218ULL);
	~ParallelMT();

	// access generators via function call with index checking
	double RandomReal(const int thread) const;
    double RandomReal(const int thread, const std::vector<double> &limits) const;

protected:

	void Cycle(unsigned long long init_key[], int key_length);
	MT19937 **generator;
	int _threads;
};

} // namespace Gaia

#endif

// BELOW IS THE ACKNOWLEDGEMENT OF USE FOR THE ORIGINAL SOURCE OF
// THE 64-BIT MERSENNE TWISTER PSEUDORANDOM NUMBER GENERATOR

/*
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/
