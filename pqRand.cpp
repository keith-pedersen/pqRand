#include "pqRand.hpp"
#include <fstream>
#include <stdexcept>
#include <assert.h>

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template<class gen64_t>
void pqr::uPRNG_64_Seeder<gen64_t>::Seed()
{
	std::stringstream ss;
	{
		std::random_device randDev;
		
		// One word for every seed, space-separated
		for(size_t i = 0; i < gen64_t::state_size; ++i)
		{
			// Unfortunately, random_device outputs 32 bits. Make a 64-bit uint.
			ss << ((uint64_t(randDev()) << 32) | uint64_t(randDev())) << " ";
		}					
		ss << gen64_t::state_size; // Terminate with state_size
	}
	Seed(ss); // Bootstrap Seed(stream)
}

////////////////////////////////////////////////////////////////////////

template<class gen64_t>
void pqr::uPRNG_64_Seeder<gen64_t>::Seed(std::string const& fileName)
{
	std::ifstream file(fileName.c_str(), std::ios::in);
	
	if(not file.is_open())
	{
		throw std::runtime_error("rpq::uPRNG_64_Seeder::Seed ... <"
			+ fileName + "> ... file not found!");
	}
	
	Seed(file); // Bootstrap Seed(stream)
	file.close();
}

////////////////////////////////////////////////////////////////////////

template<class gen64_t>
void pqr::uPRNG_64_Seeder<gen64_t>::Seed(std::istream& stream)
{
	stream >> *this; // uPRNG_64_Seeder (as a wrapper) has no state to seed
}

////////////////////////////////////////////////////////////////////////

template<class gen64_t>
void pqr::uPRNG_64_Seeder<gen64_t>::WriteState(std::string const& fileName)
{
	// CAUTION: overwrite existing without warning (ios::trunc)
	std::ofstream file(fileName.c_str(), std::ios::out | std::ios::trunc);
	
	if(not file.is_open())
	{
		throw std::runtime_error("rpq::uPRNG_64_Seeder::WriteState ... <"
			+ fileName + "> ... file cannot be created!");
	}
	
	WriteState(file); // Bootstrap WriteState(stream)
	file.close();
}

////////////////////////////////////////////////////////////////////////

template<class gen64_t>
void pqr::uPRNG_64_Seeder<gen64_t>::WriteState(std::ostream& stream)
{
	stream << *this; // uPRNG_64_Seeder (as a wrapper) has no state to write
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// BEGIN deep magic, do not touch!
uint64_t const pqr::xorshift1024_star::JUMP[pqr::xorshift1024_star::state_size] =  
{ 0x84242f96eca9c41d, 0xa3c65b8776f96855, 0x5b34a39f070b5837, 0x4489affce4f31a1e,
	0x2ffeeb0a48316f40, 0xdc2d9891fe68c022, 0x3659132bb12fea70, 0xaac17d8efa43cab8,
	0xc4cb815590989b13, 0x5ee975283d71c93b, 0x691548c86c1bd540, 0x7910c41d10a1e6a5, 
	0x0b5fc64563b3e2a8, 0x047f7684e9fc949d, 0xb99181f2d8f685ca, 0x284600e3f30e38c3};

uint64_t pqr::xorshift1024_star::operator()() 
{
	uint64_t s1;
	{
		uint64_t const s0 = state[p];
		s1 = state[p = (p + 1) & 15];
		s1 ^= s1 << 31; // a
		s1 = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
	}
	state[p] = s1;
	return s1 * UINT64_C(1181783497276652981);
}

void pqr::xorshift1024_star::Jump()
{
	uint64_t t[16] = { 0 };
	for(size_t i = 0; i < (sizeof(JUMP) / sizeof(*JUMP)); i++)
	{
		for(size_t b = 0; b < 64; b++) 
		{
			if (JUMP[i] & 1ULL << b) // What is the order of operations here?
			{
				for(size_t j = 0; j < 16; j++)
					t[j] ^= state[(j + p) & 15];
			}
			(*this)();
		}
	}

	for(size_t j = 0; j < 16; j++)
		state[(j + p) & 15] = t[j];
}
// END deep magic

////////////////////////////////////////////////////////////////////////

void pqr::xorshift1024_star::Jump(size_t nTimes)
{
	for(size_t n = 0; n < nTimes; ++n)
		Jump();
}

////////////////////////////////////////////////////////////////////////

std::ostream& pqr::operator << (std::ostream& stream, xorshift1024_star const& gen)
{
	for(size_t i = 0; i < xorshift1024_star::state_size; ++i)
		stream << gen.state[i] << " ";
	
	stream << xorshift1024_star::state_size << " ";
	
	// Store p as well (it defines the internal state too).
	stream << gen.p;
	
	return stream;
}

////////////////////////////////////////////////////////////////////////

// Seed the generator from the stream. Two formats expected (N = state_size)
// w_1 w_2 ... w_N N    --> p not specified, set to zero
// w_1 w_2 ... w_N N p  --> p specified
std::istream& pqr::operator >> (std::istream& stream, xorshift1024_star& gen)
{
	uint64_t word;
	
	for(size_t i = 0; i < xorshift1024_star::state_size; ++i)
	{	
		if(not (stream >> word))
			throw std::runtime_error("pqr::xorshift1024_star: seed stream malformed -- not enough words to fill state.");
		
		gen.state[i] = word;
	}
	
	if(not (stream >> word))
		throw std::runtime_error("pqr::xorshift1024_star: seed stream malformed -- state size not supplied.");
	else if(word not_eq xorshift1024_star::state_size)
		throw std::runtime_error("pqr::xorshift1024_star: seed stream malformed -- wrong state size.");
	
	// Read p, which exists in [0, 16). If p is not stored, then use p = 0
	if(stream >> word)
	{
		if(word > xorshift1024_star::state_size)
			throw std::runtime_error("pqr::xorshift1024_star: seed stream malformed -- p is larger than state_size");
		gen.p = word;
	}
	else gen.p = 0;
	
	return stream;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Need to instantiate the template class for the object file
template class pqr::uPRNG_64_Seeder<pqr::PRNG_t>;

////////////////////////////////////////////////////////////////////////

typename pqr::real_t pqr::pqRand_engine::RandomMantissa_Quasiuniform()
{
	uint64_t randUint = (*this)();
	
	if(randUint < minEntropy) // Add entropy to randUint until it has enough
	{
		real_t downScale = 1.; // Downscale will be set at the end
		
		{
			// We need to shift left at least once, so let's start with that
			int shiftLeft = 1; // Must use signed type, for negative exponent in pow()
			randUint <<= 1;
			
			if(randUint == 0) // Exceedingly rare, but need to handle
			{
				shiftLeft = 0; // Undo the initial left shift because ...
				
				// ... every time we draw a zero, do a 64-bit leftward shift.
				// It's like we have an infinite bit stream which we keep shifting left
				do // We already drew one zero, so we have to downscale at least once
					downScale *= scaleToU_Quasiuniform;
				while((randUint = (*this)()) == 0);
			}
			
			// Keep shifting left until the mantissa's most significant bit is
			// in the correct position
			while(randUint < minEntropy)
			{
				randUint <<= 1;
				++shiftLeft;
			}
			
			assert(shiftLeft < numBitsOfEntropyRequired);
			downScale *= std::pow(2., -shiftLeft);
						
			// Insert new bits into the gap filled by the shift left
			// Usually quite wasteful, but generally quite rare
			randUint |= ((*this)() >> (64 - shiftLeft));
		}
		
		// Set the sticky bit, to defeat round-to-even,
		// then scale by the shiftLeft, to maintain the coarse location
		return real_t(randUint | 1) * downScale;
	}
	else
		return real_t(randUint | 1);
}

////////////////////////////////////////////////////////////////////////

typename pqr::real_t pqr::pqRand_engine::RandomMantissa_Superuniform()
{
	uint64_t randUint;
	
	// Yes, this is a goto, but they are not always bad.
	// This jump spans only a small number of lines, so it's easy to read.
	generate:
		randUint = ((*this)() & bitMask_Superuniform);
		if(randUint == 0)
		{
			// Map 0 to the max mantissa, but only half the time
			// Should not cause an endless loop, unless PRNG is broken
			if(RandBool())
				randUint = maxMantissa_Superuniform;
			else
				goto generate;
		}
	
	// Equivalent code, without goto, but every call has an extra addition
	// (hence the virtue of the goto is speed).
	//~ do
		//~ randUint = (((*this)() & bitMask_Superuniform) + 1);
	//~ while((randUint == maxMantissa_Superuniform) and RandBool());
	
	return real_t(randUint);
}

////////////////////////////////////////////////////////////////////////

void pqr::pqRand_engine::Seed(std::istream& stream)
{
	// Let the PRNG seed itself, advancing the stream
	stream >> *this;
		
	// The internal state of pqRand_engine contains the bitCache,
	// which should be appended to the seed stream after PRNG details
	uint64_t word;
	if(stream >> word)
	{
		bitCache = word;
		
		// The state of bitCache can be missing, but the cacheMask cannot
		if(not (stream >> word))
			throw std::runtime_error("pqRand_engine::Seed: bitCache stored in seed, but not cacheMask");
		else
			cacheMask = word;
			
		// WARNING! What if badBits changes? What if the PRNG itself changes?
		// We assume that seeds will not be shared among different builds!
	}
	else
	{
		// If state of the bitCache is missing, simply replenish during next RandBool()
		cacheMask = replenishBitCache;
	}
}

////////////////////////////////////////////////////////////////////////

void pqr::pqRand_engine::WriteState(std::ostream& stream)
{
	// Write out the state of the generator
	stream << *this;
	
	// Now write out the state of the bitCache and the cacheMask
	stream  << " " <<  bitCache << " " << cacheMask;
}

////////////////////////////////////////////////////////////////////////

bool pqr::pqRand_engine::RandBool()
{
	if(cacheMask == replenishBitCache)
	{
		bitCache = (*this)();
		cacheMask = (uint64_t(1) << 63);
	}
	
	bool const decision = bool(cacheMask & bitCache);
	cacheMask >>= 1;
	return decision;
}

////////////////////////////////////////////////////////////////////////

void pqr::pqRand_engine::ApplyRandomSign(real_t& victim)
{
	bool const decision = RandBool();
	victim = std::copysign(victim, decision ? 1. : -1.);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqr::real_t pqr::standard_normal::operator()(pqRand_engine& gen)
{
	if(valueCached)
	{
		valueCached = false;
		return cache;
	}
	else
	{
		valueCached = true;
		two const pair = GenTwo(gen);
		cache = pair.x;
		return pair.y;
	}
}

////////////////////////////////////////////////////////////////////////

typename pqr::two pqr::standard_normal::GenTwo(pqRand_engine& gen)
{
	two pair;
	real_t u;
	
	// The Marsaglia polar method. Note that we're not doing 1 - U
	// to get a cheap random sign; this destroys the precision of U
	do
	{
		pair.x = gen.U_Q();
		pair.y = gen.U_Q();
	
		// Draw x and y from U, and reject when we don't land in the circle
		u = pair.x*pair.x + pair.y*pair.y;
		
		// Reject 2/3 of the region that rounds to 1
		// (we want the epsilon/2 to the left of 1, but not the epsilon the right of 1.
		//  There is a small region near (1, 0), (0, 1), etc. where the R region 
		//  doesn't exist, but it is vanishingly small).
		if((u == 1.) and (gen.U_Q()*3. < 2.))
			u = 2.; // Easy way to reject
	}
	while(u > 1.);

	// Give x and y a random sign via the pqRand_engine (uses its bitCache)
	gen.ApplyRandomSign(pair.x);
	gen.ApplyRandomSign(pair.y);
	
	// Implement the quantile flip-flop and scale x and y
	u = (gen.RandBool() ? 
		std::sqrt(-2.*std::log(0.5*u)/u) : 
		std::sqrt(2.*std::log1p(u/(2. - u))/u));		
	
	pair.x *= u;
	pair.y *= u;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Normal takes numbers from standard_normal and adjusts them to sigma and mu
typename pqr::two pqr::normal::GenTwo(pqRand_engine& gen)
{
	two pair = pqr::standard_normal::GenTwo(gen);
	
	pair.x = mu_ + sigma_ * pair.x;
	pair.y = mu_ + sigma_ * pair.y;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// Log_normal also takes from standard_normal, then exponentiates
typename pqr::two pqr::log_normal::GenTwo(pqRand_engine& gen)
{
	// Draw from standard normal, apply mu as a multiplicative scale
	two pair = pqr::standard_normal::GenTwo(gen);
	
	pair.x = muScale * std::exp(sigma_ * pair.x);
	pair.y = muScale * std::exp(sigma_ * pair.y);
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// When we implement the quantile flip-flop, we must draw from HalfU
typename pqr::real_t pqr::weibull::operator()(pqRand_engine& gen)
{
	real_t const hu = gen.HalfU_Q();
	
	if(gen.RandBool())
		return lambda_ * std::pow(-std::log(hu), kRecip);
	else
		return lambda_ * std::pow(std::log1p(hu/(1.-hu)), kRecip);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqr::real_t pqr::pareto::operator()(pqRand_engine& gen)
{
	return xm_ * std::pow(gen.U_Q(), negRecipAlpha);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqr::real_t pqr::exponential::operator()(pqRand_engine& gen)
{
	real_t const hu = gen.HalfU_Q();
	
	if(gen.RandBool())
		return -mu_ * std::log(hu);
	else
		return mu_ * std::log1p(hu/(1.-hu));
}
