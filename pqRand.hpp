// .o88o. .o88o. |  _ \    __ _   _ __     __| |
// 8    8 8    8 | |_) |  / _` | | '_ \   / _` |
// 8    8 8    8 |  _ <  | (_| | | | | | | (_| |
// 88oo8' `8oo88 |_| \_\  \__,_| |_| |_|  \__,_|
// 8           8 
// 8           8 
//
// pqRand: The precise quantile flip-flop random package
//   version 0.2 (Apr 2017)
// By Keith Pedersen (Keith.David.Pedersen @ gmail.com)
//   https://github.com/keith-pedersen/pqRand
//
// Special thanks to Andrew Webster, Zack Sullivan and Sebastiano Vigna.
//
// This package is designed to improve sampling from the 
// TAILS of the following distributions:
//	   * normal
//    * log-normal
//    * exponential
//    * weibull
//    * pareto
//
// This package uses C++11. Hopefully by now (2017), it is supported at your site.
// My apologies if it is not, but there are SO many reasons to use it.
//
// The principle theory behind this package is presented in the paper
// "Conditioning your quantile function," which is cited in the <copyright notice>.
// Distributed with the package are a selection of distributions which possess 
// an analytic quantile function (or, as is the case for the normal and log-normal, 
// which can be drawn using a quantile function of an underlying distribution). 
// The main methods implemented here are the quantile flip-flop and
// a "quasiuniform sample" of the semi-inclusive unit interval.
// For quality control, the user is forced to use the supplied PRNG 
// to generate uniform, unsigned integers.
//
// Development record
// =====================================================================
// version 0.1 ===> 06 Mar 17 
//		known as rTailor
// version 0.2 ===> 25 Apr 17 
//		changed name to pqRand, pushed paper to arXiv, pushed code to GitHub
//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// BEGIN <copyright notice>
//
// pqRand: The precise quantile flip-flop random package
//    Copyright(c) 2017 Keith Pedersen
//    https://github.com/keith-pedersen/pqRand
//
// Based upon work presented in:
//    "Conditioning your quantile function" <arXiv:1704.xxxxx>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistribution of source code must retain the entirety of the
//    <copyright notice>. Please cite the paper in the documentation and/or 
//    other materials provided with the distribution.
//  * Redistribution of methods in binary form (even in part) must 
//		reproduce the entire <copyright notice> in the documentation and/or 
//    other materials provided with the distribution.
//  * Neither the name of the copyright owner nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE (pqRand) IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED. ALTHOUGH THE AUTHORS MADE
//  EVERY EFFORT TO ENSURE THE FITNESS OF THIS CODE FOR ITS STATED PURPOSE, 
//  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
//  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
//  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
//  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// END <copyright notice>
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// Features to add:
//   * If a regular sample of U produces numbers with weight == 1., 
//     a sample from U that samples from the tail (u -> 0)
//     given a supplied weight << 1.

// Bugs to fix:

#ifndef PQ_RAND
#define PQ_RAND

#include <limits> // numeric_limits
#include <cmath> // exp
#include <random> // random_device / mt19937_64
#include <sstream>
#include <array>

namespace pqr
{
	// uPRNG_64_Seeder is a template class for seeding a PRNG of uint64_t.
	// It forces an initial seed AS LARGE as the state of the PRNG
	// (i.e. you can't seed 20 kB of the Mersenne twister with a 32 bit number).
	// The default seed method (no arguments) seeds from std::random_device 
	// (which internally uses /dev/urandom on many linux distros).
	// After a default seeding, one can write the generator's state to a file
	// (or a std::ostream) to allow future re-seeding from a known state.
	//
	// One can also supply a custom seed (filled in some arbitrary way)
	// by supplying a fileName/istream in the correct format
	// (see the comments above Seed/WriteSeed).
	// WARNING: The ONLY check made on the actual seed words themselves is
	// that they can actually be parsed as uint64_t, so take care 
	// NOT to supply a bunch of zeroes.
	template<class gen64_t>
	class uPRNG_64_Seeder : public gen64_t
	{
		static_assert(gen64_t::word_size == 64, "uPRNG_64_Seeder requires a 64-bit generator");
		
		protected:
			// This ctor skips the initial seed at base class level.
			// It is used by a derived class that needs to access its own 
			// virtual Seed(istream) during its initial seeding,
			// because you can't call a virtual function in a base-class ctor.
			// This is a kludge, but I don't know a better way
			// (at least not one which respects DRY [don't repeat yourself];
			//  Why do I want to keep my code DRY? It keeps the bugs away).
			explicit uPRNG_64_Seeder(bool dummy) {}
					
		public:
			// Force an initial seed, using the version of Seed taking the 
			// supplied arguments. This is an example of "perfect forwarding".
			template <typename... Args>
			uPRNG_64_Seeder(Args&&... args)
			{
				Seed(std::forward<Args>(args)...);
			}
		 
			// The format of the seed file/stream mimics std::mt19937:
			// A single line, with each word of seed space separated, 
			// terminating with the seed_size itself
			//   s_1 s_2 s_3 ... s_(state_size) state_size
			// Example for PRNG with state_size == 3
			//   8602057372317241997 13337802638347508439 10520579822156476364 3
			// Note that the seed MUST fill the entire state of the PRNG
			// ( (# of words) == seed_size),  otherwise an exception is thrown.
			
			void Seed(); // seed from random_device
			void Seed(std::string const& fileName); // seed from a file
			virtual void Seed(std::istream& stream); // seed from an input-stream
			
			void WriteState(std::string const& fileName); // write state to seed file
			virtual void WriteState(std::ostream& stream); // write state to a stream
			
			// If a derived class changes the format of the Seed stream
			// (e.g. because it adds pieces to the internal state),
			// only the virtual istream/ostream methods need redefinition. 
			// The non-virtual methods are DRY -- they use the 
			// virtual methods to do the actual seeding.
			// Note, however, that if the virtual methods are redefined, 
			// the non-virtual members must be  explicitly exposed by 
			// the derived class in order for ADL (argument dependent lookup)
			// to find them (see the "using" statements in pqRand_engine).
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	// An implementation of the xorshift1024* 
	//    By Sebastiano Vigna (vigna@acm.org), Copywrite(c) 2014
	//    http://xoroshiro.di.unimi.it/xorshift1024star.c
	// Implementation by Keith Pedersen, with API (minimally) copying std::mt19937_64
	//
	// This implementation was tested using dieharder 3.31.1, 
	// 	./xorshift2014star_tester.x | dieharder -g 200 -a -k 2
	// and passed with flying colors (no failures).
	// Sebsatiano states:
	// "... The three lowest bits of this generator are LSFRs, and thus
	// they are slightly less random than the other bits. 
	// We suggest to use a sign test to extract a random Boolean value."
	// When I tested the last 3 bits with dieharder,
	// (concatenating the last-3 from three sequential calls, x = 21321321)
	// 	./xorshift2014star_last3bits_tester.x | dieharder -g 200 -a
	// they also passed. So perhaps they're not as bad as claimed.
	// But for now, we conservatively ignore the last 3 bits in pqRand_engine
	//
	// Per Vigna: "The state must be seeded so that it is not everywhere zero."
	class xorshift1024_star
	{
		public:
			size_t static constexpr word_size = 64;
			size_t static constexpr state_size = 16;
		
		 private:
			// Lots of deep magic here, not much to comment
			uint64_t static const JUMP[state_size];
			
			// Use std::array to allow implicit assignment (operator =)
			// of this class, to quickly copy states of the generator,
			// while remaining compatible with C-style array syntax
			std::array<uint64_t, state_size> state; 
			uint64_t p;
			
			// Jump forward by 2**512 calls to the generator
			// WARNING: This seems like much deeper magic than the generator itself,
			// so I cannot guarantee that this actually works
			// (e.g. there is no dieharder test for the size of generator jumps)
			// but we can show that Jump() is commutative 
			// ((operator()(), Jump(); operator()()) = (Jump(); operator()(), operator()()))
			void Jump();
			
		public:
			// Initialize p, but not state, placing the generator in a 
			// valid, but undefined state (unless the state happens to be all zeros, 
			// in which case the generator is NOT in a valid state).
			xorshift1024_star():p(0) {}
			
			uint64_t operator()(); // generate the next number
			
			// Jump forward by (nTimes * 2**512) calls to the generator,
			// without actually calling the generator that many times.
			// Useful for parallel threads (each with its own generator)
			// so that there is no overlap in their sequences.			
			void Jump(size_t nTimes);
			
			// Declare these functions as friends
			
			// Output the state of the generator to a stream, mimicking the 
			//	stream format of std::mt19937_64 ...
			// 	s_1 s_2 s_3 ... s_(state_size) state_size p
			friend std::ostream& operator << (std::ostream& stream, xorshift1024_star const& gen);
			
			// Suck in the state of the generator from a stream, 
			// assuming the format of the above operator << .
			// Throws an exception if the stream runs out before the state is full,
			// or if the stream does not seem to follow the assumed format.
			// However, p CAN be missing (e.g. uPRNG_64_Seeder::Seed() doesn't supply p).
			// In this case, we use p = 0.
			friend std::istream& operator >> (std::istream& stream, xorshift1024_star& gen);	
	};
		
	// Actually declare the friend functions for xorshift1024_star
	// (previously declared as friends, but not as actual functions).
	std::ostream& operator << (std::ostream& stream, xorshift1024_star const& gen);
	std::istream& operator >> (std::istream& stream, xorshift1024_star& gen);	
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
			
	// mt19937 has good statistics, but its state is too large for massive multi-parallelism
	// Nonetheless, it is fully compatible with this package.
	typedef xorshift1024_star PRNG_t;
	//~ typedef std::mt19937_64 PRNG_t; 
	
	// The authors of pqRand have carefully chosen a PRNG with a 
	// sufficiently large state (but not too large), 
	// excellent uniformity properties, rigorously tested by others and 
	// vetted by the authors using the dieharder battery of tests.
	// The motivation of the pqRand package essentially requires that
	// freedom in chosing a PRNG be removed from the user, since such freedom would 
	// permit them to unwittingly destroy the benefits of pqRand.
	// Nonetheless, it is still possible to override this decision by 
	// redefining PRNG_t, but the substituted class must use the same API:
	// static fields storing the word_size and state_size, a default ctor, 
	// operator()() to return the next uint64_t, and friend operators >> and << 
	// to input/output the generator's internal state [for seeding],
	// with these operators using the same stream format as std::mt19937_64.
	// This means that std::mt19937_64 (the model for the xorshift102_star API)
	// will work as a drop-in replacement (simply switch the typedef above)
	
	// This typedef is for internal testing; there is no reason to use float/binary32
	typedef double real_t;
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
			
	// pqRand_engine wraps the chosen PRNG in an API for the quantile flip-flop,
	// providing access to:
	//   * A quasiuniform sample of U(0, 1]
	//   * A quasiuniform sample of U(0, 1/2]
	//   * A superuniform sample of U(0, 1]
	//   * A superuniform sample of U(0, 1/2]
	//   * A random bool that doesn't waste random bits
	//   * A method to efficiently apply a random sign to a floating-point
	// Future changes to the API are not anticipated, but may still be possible
	class pqRand_engine : public uPRNG_64_Seeder<PRNG_t>
	{
		private:
			// Scale the mantissa into the correct interval
			real_t static constexpr scaleToU_Quasiuniform = 
				// If conversion to real_t forces rounding, then adding 1 can't change the value
				1./(real_t(std::numeric_limits<uint64_t>::max()) + 1.);
			real_t static constexpr scaleToHalfU_Quasiuniform = 
				0.5 * scaleToU_Quasiuniform;
				
			real_t static constexpr scaleToU_Superuniform = 
				// epsilon = 2**-(P-1), we want 2**-P
				0.5 * std::numeric_limits<real_t>::epsilon();
			real_t static constexpr scaleToHalfU_Superuniform = 
				0.5 * scaleToU_Superuniform;
					
			// The last three bits of xorshift1024_star are not as good to use, per S. Vigna
			uint64_t static constexpr badBits = 3;
			uint64_t static constexpr replenishBitCache = (uint64_t(1) << (badBits - 1));
			// NOTE: if badBits == 0, replenishBitCache == 0 (a shift left of -1).
			// This may output a warning, but should still be valid.
			
			uint64_t static constexpr numBitsMantissa = std::numeric_limits<real_t>::digits;
			uint64_t static constexpr bitMask_Superuniform = 
				(uint64_t(1) << numBitsMantissa) - 1;
			uint64_t static constexpr maxMantissa_Superuniform = bitMask_Superuniform + 1;
			
			// We need to fill the mantissa, with 2 bits for rounding:
			// a buffer bit and a sticky bit. The buffer bit must use 
			// a good bit of entropy, but the sticky bit is always set to 1, 
			// so we can use the last bad bit as the sticky bit
			uint64_t static constexpr numBitsOfEntropyRequired = 
				numBitsMantissa + 1 + ((badBits > 0) ? (badBits - 1): 1);
			// When the random uint is less than minEntropy, we need more entropy.
			// minEntropy need only skip (badBits - 1), 
			uint64_t static constexpr minEntropy = 
				(uint64_t(1) << (numBitsOfEntropyRequired - 1));
						
			uint64_t bitCache; // A cache of random bits for RandBool
			uint64_t cacheMask; // Selects one bit from bitCache
		
			// Draw a random, quasiuniform mantissa from (0, 2**53],
			// with 2**53 being half as probable as 2**53-1
			real_t RandomMantissa_Quasiuniform();
			// Draw a random, superuniform mantissa from (0, 2**53],
			// with 2**53 being half as probable as 2**53-1
			real_t RandomMantissa_Superuniform();
			
		public:
			// Tell the base class not to seed, because we need to call our 
			// newly redefined virtual Seed() for the intiial seed.
			template <typename... Args>
			pqRand_engine(Args&&... args):uPRNG_64_Seeder(false)
			{
				Seed(std::forward<Args>(args)...);
			}			 
			
			// Redefine the base class virtuals, because we need to 
			// store/refresh the state of the bitCache when we write/seed
			// Also expose the base class functions of the same name, 
			// otherwise they will be hidden during ADL.
			virtual void Seed(std::istream& stream);
			using uPRNG_64_Seeder<PRNG_t>::Seed;
			
			virtual void WriteState(std::ostream& stream);
			using uPRNG_64_Seeder<PRNG_t>::WriteState;
					
			bool RandBool(); // An ideal coin flip
			void ApplyRandomSign(real_t& victim); // Give victim a random sign (+/-)
			
			// Access the semi-inclusive or inclusive unit inveral, 
			// with inclusive endpoints being half as probable as the their neighbors
			// (for use in the quantile flip-flop)
			// U => 		U(0, 1]
			// HalfU => U(0, 1/2]
				
			inline real_t U_Q() 
			{
				return scaleToU_Quasiuniform * RandomMantissa_Quasiuniform();
			}
			
			inline real_t HalfU_Q() 
			{
				return scaleToHalfU_Quasiuniform * RandomMantissa_Quasiuniform();
			}
			
			inline real_t U_S() 
			{
				return scaleToU_Superuniform * RandomMantissa_Superuniform();
			}
			
			inline real_t HalfU_S()
			{
				return scaleToHalfU_Superuniform * RandomMantissa_Superuniform();
			}
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
		
	// With the worker classes defined, we can finally use the
	// quantile flip-flop to define some distributions.
	// Luckily, we've already done all the hard work!
	
	// A simple struct to store a pair of objects, used by the 
	// normal distributions, since Marsaglia polar generates two number per call
	struct two
	{
		real_t x;
		real_t y;
		
		two(real_t const x_in, real_t const y_in):
			x(x_in), y(y_in) {}
		
		two() {} // No default initialization, be careful
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	// Sample from the standard normal distribution using the Marsaglia polar method 
	// (this implementation is modeled after GNU's std::normal_distribution,
	// but with the added precision of the quantile flip-flop and U_Q).
	class standard_normal
	{
		private:
			real_t cache; // Marsaglia polar samples two, cache one if we only request one
			bool valueCached; // if(valueCached == true), cache has the next value
			
		public:
			standard_normal():
				valueCached(false) {}
			
			real_t operator()(pqRand_engine& gen); // Return a single value
			virtual two GenTwo(pqRand_engine& gen); // Return a pair of values
			// GenTwo is virtual because we can re-use the caching functionality
			// for normal and log_normal by simply redefining GenTwo
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	// Bootstrap normal via standard_normal ... x = mu + sigma * x_standard
	class normal : public standard_normal
	{
		protected:
			real_t const mu_;
			real_t const sigma_;
			
		public:
			normal(real_t const mu_in, real_t const sigma_in):
				standard_normal(), mu_(mu_in), sigma_(sigma_in) {}
			
			// Overwrite GenTwo, to add sigma and mu
			virtual two GenTwo(pqRand_engine& gen);
			
			inline real_t mu() {return mu_;}
			inline real_t sigma() {return sigma_;}
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	// Bootstrap log_normal from normal ... x = e**mu * e**(sigma*x_standard)
	class log_normal : public normal
	{
		private:
			real_t const muScale; // Exponentiate mu once, multiply times all returned values
			// This is done to prevent assumed addition cancellation in 
			// exp(mu + sigma * x), but is the cancellation actually avoided?
			
		public:
			log_normal(real_t const mu_in, real_t const sigma_in):
				normal(mu_in, sigma_in),
				muScale(std::exp(mu_)) {}
				
			virtual two GenTwo(pqRand_engine& gen);
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	// Weibull has an invertible CDF: x = lambda*(-log(1-u))**(1/k)
	// Use a quantile flip-flop and U_Q
	class weibull
	{
		private:
			real_t const lambda_;
			real_t const k_;
			real_t const kRecip;
			
		public:
			weibull(real_t const lambda_in, real_t const k_in):
				lambda_(lambda_in), k_(k_in), kRecip(1./k_) {}
				
			real_t operator()(pqRand_engine& gen);
			
			inline real_t lambda() {return lambda_;}
			inline real_t k() {return k_;}
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	// Pareto has an invertible CDF: x = x_m * (1 - u)**(-1/alpha)
	// Doesn't need a flip-flop if we remove floating-point cancellation and use U_Q
	class pareto
	{
		private:
			real_t const xm_;
			real_t const alpha_;
			real_t const negRecipAlpha;
			
		public:
			pareto(real_t const xm_in, real_t const alpha_in):
				xm_(xm_in), alpha_(alpha_in), negRecipAlpha(-1./alpha_) {}
				
			real_t operator()(pqRand_engine& gen);
			
			inline real_t xm() {return xm_;}
			inline real_t alpha() {return alpha_;}
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	// exponential has an invertible CDF: x = -mu * log(1 - u)
	class exponential
	{
		private:
			real_t const mu_;
			
		public:
			exponential(real_t const mu_in):
				mu_(mu_in) {}
				
			real_t operator()(pqRand_engine& gen);
			
			inline real_t mu() {return mu_;}
	};	
};

#endif
