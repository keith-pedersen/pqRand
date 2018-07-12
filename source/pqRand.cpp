/* pqRand: The precise quantile random package
 * Copyright (C) 2017 Keith Pedersen (Keith.David.Pedersen@gmail.com)
 * 
 * This package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This package is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the COPYRIGHT_NOTICE for more details.
 * 
 * Under Section 7 of GPL version 3, you are granted additional
 * permissions described in the GCC Runtime Library Exception, version
 * 3.1, as published by the Free Software Foundation.
 * 
 * You should have received a copy of the GNU General Public License and
 * a copy of the GCC Runtime Library Exception along with this package;
 * see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
 * <http://www.gnu.org/licenses/>.
*/

#include "../include/pqRand.hpp"
#include <fstream>
#include <assert.h>
#include <sstream>
#include <random> // random_device, mt19937

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template<class prng_t>
void pqRand::seeded_uPRNG<prng_t>::Seed()
{
	// Construct the minimal state-string into a temporary stream
	std::stringstream ss;
	{
		std::random_device randDev;
		
		// Insert one word for every state_size, space-separated
		switch(seeded_uPRNG::word_size)
		{
			case 32:
				for(size_t i = 0; i < prng_t::state_size; ++i)
				{
					ss << randDev() << " ";
				}
			break;
			
			case 64:
				for(size_t i = 0; i < prng_t::state_size; ++i)
				{
					// Unfortunately, random_device outputs 32 bits. Make a 64-bit uint.
					ss << ((uint64_t(randDev()) << 32) bitor uint64_t(randDev())) << " ";
				}
			break;
		}
	}
	
	ss << prng_t::state_size; // Terminate with state_size
	this->Seed_FromStream(ss);
}

////////////////////////////////////////////////////////////////////////

template<class prng_t>
void pqRand::seeded_uPRNG<prng_t>::Seed_FromFile(std::string const& filePath)
{
	std::ifstream file(filePath.c_str(), std::ios::in);
	
	if(not file.is_open())
	{
		throw std::ifstream::failure("pqRand::seeded_uPRNG::Seed ... seed file <"
			+ filePath + "> ... cannot be opened (probably does not exist)!");
	}
	
	this->Seed_FromStream(file);
	file.close();
}

////////////////////////////////////////////////////////////////////////

template<class prng_t>
void pqRand::seeded_uPRNG<prng_t>::Seed_Reuse(std::string const& filePath)
{
	std::ifstream file(filePath.c_str(), std::ios::in);
	
	if(not file.is_open())
	{
		Seed();
		WriteState(filePath);
	}
	else
	{
		this->Seed_FromStream(file);
		file.close();
	}
}

////////////////////////////////////////////////////////////////////////

template<class prng_t>
void pqRand::seeded_uPRNG<prng_t>::Seed_FromString(std::string const& seedString)
{
	std::stringstream stream(seedString);
	this->Seed_FromStream(stream);
}

////////////////////////////////////////////////////////////////////////

template<class prng_t>
void pqRand::seeded_uPRNG<prng_t>::Seed_FromStream(std::istream& stream)
{
	// Nothing to do but put the stream into the PRNG, because 
	// seeded_uPRNG (as a wrapper) has no other state to seed
	stream >> *this; 
}

////////////////////////////////////////////////////////////////////////

template<class prng_t>
void pqRand::seeded_uPRNG<prng_t>::WriteState(std::string const& filePath)
{
	// CAUTION: overwrite existing file without warning (ios::trunc)
	std::ofstream file(filePath.c_str(), std::ios::out | std::ios::trunc);
	
	if(not file.is_open())
	{
		throw std::ifstream::failure("pqRand::seeded_uPRNG::WriteState ... seed file <"
			+ filePath + "> ... cannot be created or overwritten!");
	}
	
	this->WriteState_ToStream(file);
	file.close();
}

////////////////////////////////////////////////////////////////////////

template<class prng_t>
std::string pqRand::seeded_uPRNG<prng_t>::GetState()
{
	std::stringstream string;
	this->WriteState_ToStream(string);
	return string.str();
}

////////////////////////////////////////////////////////////////////////

template<class prng_t>
void pqRand::seeded_uPRNG<prng_t>::WriteState_ToStream(std::ostream& stream)
{
	// Nothing to do but put the stream into the PRNG, because 
	// seeded_uPRNG (as a wrapper) has no other state to add 
	stream << *this;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// BEGIN deep magic, do not touch!
// xorshift1024_starPhi (SV's 2017 edit)
uint64_t pqRand::xorshift1024_star::operator()() 
{
	uint64_t s1;
	{
		uint64_t const s0 = state[p];
		s1 = state[p = (p + 1) & 15]; // Fast modulo 16 (that's an & not a %)
		s1 ^= s1 << 31; // a
		s1 = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
	}
	state[p] = s1;
	return s1 * 0x9e3779b97f4a7c13;
}

void pqRand::xorshift1024_star::Jump()
{
	static constexpr uint64_t JUMP[] = { 0x84242f96eca9c41d,
		0xa3c65b8776f96855, 0x5b34a39f070b5837, 0x4489affce4f31a1e,
		0x2ffeeb0a48316f40, 0xdc2d9891fe68c022, 0x3659132bb12fea70,
		0xaac17d8efa43cab8, 0xc4cb815590989b13, 0x5ee975283d71c93b,
		0x691548c86c1bd540, 0x7910c41d10a1e6a5, 0x0b5fc64563b3e2a8,
		0x047f7684e9fc949d, 0xb99181f2d8f685ca, 0x284600e3f30e38c3
	};
	
	uint64_t t[16] = { 0 };
	for(uint64_t i = 0; i < (sizeof(JUMP) / sizeof(*JUMP)); i++)
	{
		for(uint64_t b = 0; b < 64lu; b++) 
		{
			if (JUMP[i] & UINT64_C(1) << b) // What is the order of operations here?
			{
				for(uint64_t j = 0; j < 16lu; j++)
					t[j] ^= state[(j + p) & 15lu];
			}
			(*this)();
		}
	}

	for(uint64_t j = 0; j < 16lu; j++)
		state[(j + p) & 15lu] = t[j];
}
// END deep magic

////////////////////////////////////////////////////////////////////////

// Write the state to the stream
std::ostream& pqRand::operator << (std::ostream& stream, xorshift1024_star const& gen)
{
	// Put each word in the stream
	for(size_t i = 0; i < xorshift1024_star::state_size; ++i)
		stream << gen.state[i] << " ";
	
	// Add the state_size
	stream << xorshift1024_star::state_size << " ";
	
	// Store p as well (it defines the internal state too).
	stream << gen.p;
	
	return stream;
}

////////////////////////////////////////////////////////////////////////

// Seed the generator from the stream. Two formats expected (N = state_size)
// s_1 s_2 ... s_N  N    --> p not specified, set to zero
// s_1 s_2 ... s_N  N p  --> p specified
std::istream& pqRand::operator >> (std::istream& stream, xorshift1024_star& gen)
{
	uint64_t word;
	
	for(size_t i = 0; i < xorshift1024_star::state_size; ++i)
	{	
		if(not (stream >> word))
			throw pqRand::seed_error("pqRand::xorshift1024_star: seed stream malformed -- not enough words to fill state.");
		
		gen.state[i] = word;
	}
	
	if(not (stream >> word))
		throw pqRand::seed_error("pqRand::xorshift1024_star: seed stream malformed -- state size not supplied.");
	else if(word not_eq xorshift1024_star::state_size)
		throw pqRand::seed_error("pqRand::xorshift1024_star: seed stream malformed -- wrong state size.");
	
	// Read p, which exists in [0, 16). If p is not stored, then use p = 0
	if(stream >> word) // If there is another parse-able word, it must be p
	{
		if(word >= xorshift1024_star::state_size)
			throw pqRand::seed_error("pqRand::xorshift1024_star: seed stream malformed -- p is larger than state_size");
		gen.p = word;
	}
	else gen.p = 0;
	
	return stream; // There might be more state to read
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Need to instantiate the template class for the object file or shared library
template class pqRand::seeded_uPRNG<pqRand::PRNG_t>;

////////////////////////////////////////////////////////////////////////

// When randUint does not have enough entropy, we make sure it has P+2 bits
typename pqRand::real_t pqRand::engine::U_uneven_TopUpEntropy(result_type randUint)
{
	// downScale reverses the leftward shift, so the uniform variate doesn't move
	// We need to shift randUint left at least once, so we start with that
	
	real_t downScale = real_t(0.5) * scaleToU_uneven;
	{
		size_t shiftLeft = 1; // Must use signed type, for negative exponent in exp2()
		randUint <<= 1;
		
		if(randUint == 0) // Exceedingly rare, but need to handle
		{
			shiftLeft = 0; // Undo the initial left shift because ...
			downScale = scaleToU_uneven;
			
			// ... every time we draw a zero, do a 64-bit leftward shift.
			// It's like we have an infinite bit stream which we keep shifting left
			do // We already drew one zero, so we have to downscale at least once
				downScale *= scaleToU_uneven;
			while((randUint = (*this)()) == 0);
		}
		
		// Keep shifting left until the mantissa's most significant bit is
		// in the correct position
		while(randUint < minEntropy)
		{
			randUint <<= 1;
			++shiftLeft;
			downScale *= real_t(0.5);
		}
		
		// This has been tested enough
		//~ assert(size_t(shiftLeft) < numBitsOfEntropyRequired);
		
		// Using downcale is empirically faster than exp2(-shiftLeft)
		// WARNING: this requires shiftLeft to be an int, so it can be negative
		//~ downScale = std::exp2(-shiftLeft); 
					
		// Insert new bits into the gap filled by the shift left
		// Usually quite wasteful, but generally rare enough
		randUint or_eq ((*this)() >> (numBitsPRNG - shiftLeft));
	}
	
	// Make randUint odd, to defeat round-to-even,
	// then downScale, to maintain the coarse location
	return real_t(randUint bitor result_type(1)) * downScale;
}

////////////////////////////////////////////////////////////////////////

void pqRand::engine::Seed_FromStream(std::istream& stream)
{
	// Seed the base class, advancing the stream
	seeded_uPRNG<PRNG_t>::Seed_FromStream(stream);
		
	// The internal state of pqRand_engine contains the bitCache,
	// which should be appended to the seed stream after PRNG details
	result_type word;
	if(stream >> word)
	{
		bitCache = word;
		
		// The bitCache can be missing; but if the bitCache is there, 
		// the cacheMask cannot be missing
		if(not (stream >> word))
			throw pqRand::seed_error("pqRand::engine::Seed: bitCache stored in seed, but not cacheMask");
		else
			cacheMask = word;
			
		// WARNING! What if badBits changes? What if the PRNG itself changes?
		// We assume that seeds will not be shared among different builds!
	}
	else
		DefaultInitializeBitCache();
}

////////////////////////////////////////////////////////////////////////

void pqRand::engine::WriteState_ToStream(std::ostream& stream)
{
	// Write out the state of the underlying PRNG
	seeded_uPRNG<PRNG_t>::WriteState_ToStream(stream);
	
	// Now write out the state of the bitCache and the cacheMask
	stream  << " " <<  bitCache << " " << cacheMask;
}

////////////////////////////////////////////////////////////////////////

void pqRand::engine::DefaultInitializeBitCache()
{
	// By setting cacheMask to this value, we ensure that the next call to RandBool()
	// will induce the bitCache to be replenished and reset.
	cacheMask = replenishBitCache;
	// Even though we will replenish the bit cache the first time access is attempted, 
	// null initialize the bit cache to prevent spurious valgrind complaints.
	bitCache = 0;
}

////////////////////////////////////////////////////////////////////////

bool pqRand::engine::RandBool()
{
	// The cacheMask starts at the leftmost bit and moves right
	if(cacheMask == replenishBitCache)
	{
		// When the cacheMask has moved too far right ... 
		bitCache = (*this)(); // Get a new set of random bits
		cacheMask = (result_type(1) << (numBitsPRNG - 1)); // Reset the cacheMask
	}
	
	bool const decision = bool(cacheMask bitand bitCache);
	cacheMask >>= 1;
	return decision;
}
