// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Library/Random.cc
//
// This source file contains the definitions for the MT19937 object and the
// wrapper class ParallelMT.

#include <Random.hpp>
#include <Exception.hpp>

namespace Gaia {

// construct via single seed
MT19937::MT19937 ( unsigned long long seed ) {

    mti = NN + 1;

    mt[0] = seed;

    for (mti = 1; mti < NN; mti++)
        mt[mti] = (6364136223846793005ULL *
			(mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

// construct via array
MT19937::MT19937(unsigned long long init_key[], unsigned long long key_length){

    unsigned long long i, j, k;

    mt[0] = 19650218ULL;

    for (mti = 1; mti < NN; mti++) {

        mt[mti] = (6364136223846793005ULL *
			(mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
    }

    i = 1; j = 0;

    k = (NN>key_length ? NN : key_length);

    for (; k; k--) {

        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
        + init_key[j] + j; /* non linear */
        i++; j++;

        if (i >= NN) {

            mt[0] = mt[NN-1]; i=1;
        }

        if (j >= key_length) {

            j=0;
    	}
	}

    for (k = NN - 1; k; k--) {

		// non linear
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) *
			2862933555777941757ULL)) - i;

        i++;

        if (i >= NN) { mt[0] = mt[NN-1]; i = 1; }
    }

    mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
}

// generates a random integer on [0, 2^64-1]-interval
unsigned long long MT19937::RandomInteger( void ) {

	int i;
	unsigned long long x;
	unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}

// generates a random number on [0,1]-real-interval
double MT19937::RandomReal() {
    return (RandomInteger() >> 11) * (1.0/9007199254740991.0);
}

// construct MT19937 family from first seed
ParallelMT::ParallelMT(const int threads, const unsigned long long first_seed){

	// set the `generator` pointer to nullptr immediately
	generator = nullptr;

	// save `threads` for destructor
	_threads = threads;

	// parent generator
	MT19937 parent(first_seed);

	// seed array has length `threads` for derived generators
	unsigned long long init_key[threads];
	for (int i = 0; i < threads; i++)
		init_key[i] = parent.RandomInteger();

	// build set of parallel generators
	generator = new MT19937 *[threads];
	for (int i = 0; i < threads; i++){

		// cycle the seed array
		Cycle(init_key, threads);

		// initialize the generator
		generator[i] = new MT19937(init_key, threads);
	}
}

// modify the seed array for the next generator
void ParallelMT::Cycle(unsigned long long init_key[], int key_length){

	for ( int i = 0; i < key_length; i++){

		if (init_key[i] < 18446744073709551615ULL - 2){
			// adding a small number results in maximal change
			init_key[i] += 2;

		} else {
			// cycle value around if beyond expected limit
			init_key[i] -= 18446744073709551615ULL - 2;
		}
	}
}

// delete the generators
ParallelMT::~ParallelMT(){

	if (generator){

		for (int i = 0; i < _threads; i++)
			delete generator[i];

		delete[] generator;
		generator = nullptr;
	}
}

// accessor for generators
double ParallelMT::RandomReal(const int thread) const {

	if (thread < _threads)
		return generator[thread] -> RandomReal();

	else throw IndexError("From ParallelMT::RandomReal(), you requested "
    "a generator by thread number that was out of bounds!");
}

// accessor for generators
double ParallelMT::RandomReal(const int thread,
    const std::vector<double> &limits) const {

    if ( limits.size() != 2 ) throw IndexError("From ParallelMT::RandomReal(), "
    "the `limits` vector size must exactly equal 2!");

	if (thread < _threads) return limits[0] + ( (limits[1] - limits[0]) *
        generator[thread] -> RandomReal() );

	else throw IndexError("From ParallelMT::RandomReal(), you requested "
    "a generator by thread number that was out of bounds!");
}

} // namespace Gaia
