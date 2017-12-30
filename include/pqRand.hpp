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

/*! @file pqRand.hpp
 *  @brief Establishes pqRand and defines the PRNG engine needed by distributions
 * 
 *  The important methods implemented here are the PRNG seeding class (pqRand::seeded_uPRNG),
 *  the underlying PRNG (pqRand::xorshift1024_star) its the pqRand::engine which wraps them together.
 * 
 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
 *  @date 2017
*/

#ifndef PQ_RAND
#define PQ_RAND

#include <limits> // numeric_limits
#include <array>
#include <vector>
#include <random> // mt19937

namespace pqRand //! @brief The namespace of the pqRand package
{	
	// Declare the @mainpage inside the namespace, so we can use unqualified names
	// (e.g. xorshift1024_star instead of pqRand::xorshift1024_star)
	
	/*! @mainpage 
	 * 
	 *  @brief pqRand
	 *  ==================================================================== 
	 *  The precise quantile random package \n
	 *  	https://github.com/keith-pedersen/pqRand \n
	 *  
	 *  \verbatim
							 ____                        _
		 .o88o. .o88o. |  _ \    __ _   _ __     __| |
		 8    8 8    8 | |_) |  / _` | | '_ \   / _` |
		 8    8 8    8 |  _ <  | (_| | | | | | | (_| |
		 88oo8' `8oo88 |_| \_\  \__,_| |_| |_|  \__,_|
		 8           8
		 8           8

		 \endverbatim
	 * 
	 *  @version 0.5.0 (beta); Dec 2017
	 *  @author  Keith Pedersen (Keith.David.Pedersen@gmail.com) \n
	 *  
	 *  @copyright 2017 (see COPYRIGHT_NOTICE)
	 * 
	 *  Special thanks to Zack Sullivan, Andrew Webster and Sebastiano Vigna.
	 * 
	 * 
	 *  Description
	 *  --------------------------------------------------------------------
	 *  
	 *  This package is designed to provide high-precision samples from the following distributions:
	 *  	+ \ref uniform
	 *  	+ \ref normal
	 *  	+ \ref log_normal "log-normal"
	 *  	+ \ref exponential
	 *  	+ \ref weibull
	 *  	+ \ref pareto
	 *  	+ \ref logistic
	 *  	+ \ref log_logistic "log-logistic"
	 *  	+ \ref gammaDist "gamma"
	 *    
	 *  Theory @anchor theory
	 *  -------
	 *  In contrast to many existing implementations of random sampling
	 *  (e.g. C++11's <\b random> suite and \b numpy.random), pqRand can deliver
	 *  the <i> most precise </i> sample possible with floating point arithmetic.
	 *  This is done via the improved inversion method described in
	 *  "Reconditioning your quantile function," https://arxiv.org/abs/1704.07949 (also distributed in doc/).
	 *  This improved method pairs an uneven uniform variate with the "quantile flip-flop".
	 *  With parallelization in mind, pqRand's default choice of
	 *  pseudo-random number generator (PRNG) is \ref xorshift1024_star,
	 *  which balances a long period with a compact state (so that each thread can have its own PRNG).
	 *  Nonetheless, the PRNG can be rather easily changed by redefining \ref pqRand::PRNG_t.
	 * 
	 *  Many of the distributions distributed with pqRand use the improved inversion method 
	 *  (i.e. they have an analytic quantile function). 
	 *  Others, like \ref normal and \ref gammaDist, use rejection sampling
	 *  The latter demonstrates that pqRand can generate a 
	 *  maximally precise sample of \em any distribution
	 *  (given a sufficiently efficient rejection algorithm).
	 * 
	 *  Usage
	 *  --------------------------------------------------------------------
	 *  1. Create one \ref engine (per thread) to draw random numbers.
	 *     - "Auto-seed" is highly recommended (see engine::Seed).
	 *     - The seed can be saved via engine::WriteState be for future audits.
	 *  2. Create as many distributions as you want.
	 *  3. Call the distributions like functions, passing the engine as the argument
	 * 	 (this mimics the API of the C++11 <\b random> suite).
	 * 
		\code
			pqRand::engine gen;  //........................... Create the random engine with auto-seed (default constructor)
			gen.WriteState("./mySeed.dat");  //............... Save the seed for future auditing
			pqRand::exponential expo(2);  //.................. Instantiate the distribution
			
			auto const one = expo(gen);  //................... Sample one variate 
			auto const many = expo.GetSample(1e6, gen);  //... Sample many variates
			
		\endcode
	 *  For more examples, see "examples/".
	 * 
	 * 
	 *  Python (pYqRand) @anchor pYqRand
	 *  -------------------------------------------------------------------
	 *  \ref pYqRand is the Python version of pqRand. The API is 99% identical to C++, 
	 *  but is independently documented via docstrings (e.g. Python "help")
	 *  \code
			>>> import pYqRand as pqr
			>>> help(pqr)
			>>> help(pqr.engine.GetState)
			>>>
			>>> #//........................................... Replicate the sample produced by the Usage snippet
			>>> gen = pqr.engine()
			>>> gen.Seed_FromFile("./mySeed.dat")  #//........ Re-seed from the stored state
			>>> expo = pqr.exponential(2.)
			>>> 
			>>> one = expo(gen)
			>>> many = expo.GetSample(1e6, gen)
		 \endcode
	 *  The minor differences (due to limitations of the Python language) 
	 *  are enumerated at the top of \ref pYqRand's docstring.	 * 
	 * 
	 *  Features to consider in future versions
	 *  --------------------------------------------------------------------
	 *  - If a regular sample of \f$ U \f$ produces numbers with \f$ weight = 1 \f$, 
	 *    supply a sample of \f$ U \f$ from the tail \f$(u \to 0)\f$,
	 *    given a supplied \f$ weight \ll 1 \f$.
	*/
	 
	class xorshift1024_star; // Forward declare this generator for the typedef below

	/*! @brief \ref PRNG_t is the PRNG used by \ref engine
	 * 
	 *  The existential purpose of the pqRand package requires
	 *  denying the user the freedom to chose the PRNG at runtime, 
	 *  since a high-precision sample from a crappy PRNG is not possible.
	 *  Nonetheless, one can change the hardcoded PRNG by 
	 *  redefining \ref PRNG_t and rebuilding; provided the replacement satisfies the 
	 *  \ref prng_requirements "\c prng_t requirements" of seeded_uPRNG.
	 *  One such candidate is \c std::mt19937_64, which was used as 
	 *  the model for the seeded_uPRNG API.
	 * 
	 *  \warning If PRNG_t is redefined, one may have to change \ref PRNG_CAN_JUMP
	 * 
	 *  The authors of pqRand have carefully chosen xorshift1024_star as the default PRNG,
	 *  since it has a "large but not too large" period (\f$ 2^{1024} \f$), 
	 *  excellent uniformity properties, has been rigorously tested by others and 
	 *  vetted by the authors using the dieharder battery of tests.
	 *  While \c std::mt19937_64 has equivalent statistics, 
	 *  its state is too large for massive multi-parallelism.
	*/
	
	/////////////////////////////////////////////////////////////////////
		
	typedef xorshift1024_star PRNG_t;
	#define PRNG_CAN_JUMP 1 //!< Does the PRNG have a Jump() function? See \ref jumping_PRNG
	
	/////////////////////////////////////////////////////////////////////
	
	//~ typedef std::mt19937_64 PRNG_t;
	//~ #define PRNG_CAN_JUMP 0 //!< Does the PRNG have a Jump() function? See \ref jumping_PRNG
	
	/////////////////////////////////////////////////////////////////////
	
	//~ typedef std::mt19937 PRNG_t;
	//~ #define PRNG_CAN_JUMP 0 //!< Does the PRNG have a Jump() function? See \ref jumping_PRNG
	
	/////////////////////////////////////////////////////////////////////
			
	//! @brief Typedef the real type to allow experimentation with float/binary32.
	//~ typedef float real_t;
	typedef double real_t;
	
	/////////////////////////////////////////////////////////////////////

	//! @brief An exception thrown when parsing state-strings fails
	class seed_error : public std::runtime_error
	{
		public:
			seed_error(std::string const& what_in):runtime_error(what_in) {}
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	inline real_t Squared(real_t const x) {return x*x;}
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
		
	/*! @brief A wrapper for a PRNG (32 or 64 bit), providing a seeding interface.
	 * 
	 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
	 * 
	 *  This class is designed to simplify seeding for \ref xorshift1024_star.
	 * 
	 *  Seeding fills the \em entire state of the PRNG
	 *  (because you can't properly seed a PRNG with 128 bytes of state using a 4-byte int).
	 *  The recommended method is Seed(), which "auto-seeds" the PRNG.
	 *  This auto-seed can then be stored for later re-seeding.
	 * 
	 *  The state of the PRNG is converted to an ASCII "state-string" by GetState(), 
	 *  or it can be written directly to a file by WriteState(). 
	 *  This permits re-seeding by Seed_FromString() or Seed_FromFile() (respectively). 
	 *  The PRNG can also be seeded from a user-generated state-string, 
	 *  provided that it is in the correct format (see Seed()).
	 *  
	 *  \anchor jumping_PRNG
	 *  \note The pre-processor macro \ref PRNG_CAN_JUMP should indicate
	 *  if PRNG_t::Jump() exists; a function which jumps the PRNG forward by a 
	 *  sufficiently huge number of calls (e.g. \f$ O(\sqrt{\text{PRNG period}}) \f$).
	 *  A Jump() function is useful for repeatably launching parallel threads, 
	 *  each with their own PRNG, but where the PRNGs sequences are 
	 *  guaranteed to be non-overlapping (and therefore not correlated). 
	 *  The two alternatives don't work:
	 *  - Feeding each thread its own random seed might might 
	 *  accidentally cause sequence collisions between a pair of threads. 
	 *   Combinatorically, this is too probable and too calamitous to ignore.
	 *  - Parallel threads can't use the same generator because 
	 *  thread scheduling is non-deterministic (and thus not repeatable).
	 * 
	 *  \warning Seed_FromFile() and Seed_FromString() do not check the format of
	 *  state-strings before forwarding them to the PRNG. If they are malformed, 
	 *  an exception may be thrown by the <c> friend operator>>(std::istream&, PRNG_t&)</c>.
	 *  
	 *  \anchor SeedWrite_Virtual
	 *  Seed_FromStream() and WriteState_ToStream() are
	 *  simple wrappers to \c prng_t 's \c operator<< and \c operator>>.
	 *  These are the actual worker functions called by the public Seed/State functions, 
	 *  and are \b virtual so that a derived PRNG class may add 
	 *  \em more information to its internal state (beyond the state of the underlying PRNG)
	 *  without redefining the Seed/State functions.
	 * 
	 *  \internal There are some extra tabs in the prng_requirement to get the formatting right in Doxygen. 
	 *            They seem to be needed for the bullets that contain the code examples. \endinternal
	 *  \anchor prng_requirements
	 * 
	 *  \param prng_t    
	 *  A PRNG class supplying \c result_type. It must:
	 *  	- have a nullary constructor 
	 *    	- possess the following static fields and functions
	 * 	\code
			 typedef result_type; // The unsigned integer type returned by the generator
			 static const word_size; // Number of bits per word (a single variate, stored in result_type)
			 static const state_size; // Number of words in the PRNG state (rounded up)	       
			 static constexpr result_type min(); // The minimum word returned by the PRNG
			 static constexpr result_type max(); // The maximum word returned by the PRNG
			\endcode
	 *    	- implement the following functions to 
	 *      input/output the state of the generator in human-readable, ASCII format
	 *  \code
			 friend operator>>(std::istream&, prng_t&); // Read the stream into the PRNG's state
			 friend operator<<(std::ostream&, prng_t&); // Write the PRNG's state into the stream
		 \endcode
	 *		- The PRNG \em must be able to repeatably configure itself from 
	 *      the "minimal" state-string described by Seed().
	*/
	template<class prng_t>
	class seeded_uPRNG : public prng_t
	{
		public: 
			// We can not document a "using" name and expect Doxygen to pick it up
			using prng_t::word_size; // Number of bits per PRNG word
			using prng_t::state_size; // Number of words in the PRNG state (rounded up).
			//using typename prng_t::result_type; // The unsigned integer type returned by the generator
						
			// Define these using lower-case syntax for comparability with std::generate_canonical
			using prng_t::min; // The smallest value this PRNG can return
			using prng_t::max; // The largest value this PRNG can return
			
			static_assert(((word_size == 32) or (word_size == 64)), 
				"pqRand::seeded_uPRNG requires either a 32-bit or 64-bit prng_t");
			
		protected:
			/*! @brief Set the state of the PRNG from a state-string stored in a stream, 
			 *  using \c prng_t's <c> friend operator>> </c>.
			 * 
			 *  See \ref SeedWrite_Virtual "the detailed class description".
			 * 
			 *  @param stream 
			 *  The state input.
			*/
			virtual void Seed_FromStream(std::istream& stream);
			
			/*! @brief Write the state-string of the PRNG to \a stream, using \c prng_t's <c> friend operator<< </c>.
			 * 
			 *  See \ref SeedWrite_Virtual "the detailed class description".
			 * 
			 *  @param stream 
			 *  The state output.
			*/
			virtual void WriteState_ToStream(std::ostream& stream);
								
		public:
			/*! @brief Construct the PRNG using its default constructor. Auto-seed if requested.
			 * 
			 *  Note that the default constructor auto-seeds. 
			 *  This scheme allows the best of both worlds:
			 *  	- Automatically use high-quality seed in most cases (lest you forget).
			 *  	- Allow seeding from stored states when explicitly requested 
			 * 		by passing false (see engine::engine for an example).
			 * 
			 *  \param autoSeed 
			 *  Perform an autoSeed (\c true) by calling Seed()
			 *  or defer seeding till later (\c false).
			 * 
			 *  \warning If seeding is deferred, the state is not initialized. Be careful.
			*/
			explicit seeded_uPRNG(bool const autoSeed = true):
				prng_t()
			{
				if(autoSeed) Seed(); // else defer seed to user
			}
			// NOTE: Access to non-auto Seed from inside the ctor is difficult.
			// states and files are both passed as strings, so which do you want?
			// Plus, Python can't overload functions, so we need a uniform ctor.
			
			virtual ~seeded_uPRNG() {}
			
			/*! @brief Auto-seed the generator using \c std::random_device.
			 * 
			 *  <a href="http://en.cppreference.com/w/cpp/numeric/random/random_device">
			 *  \c std::random_device </a> 
			 *  is a random number generator that should produce non-deterministic random numbers.
			 *  On many GNU + Linux systems, \c std::random_device uses \c /dev/urandom,
			 *  which itself uses environmental noise to seed a cryptographic PRNG. 
			 *  Seed() calls \c std::random_device until the generator's state is filled.
			 *  
			 *  <b>Minimal state-string</b>
			 * 
			 *  Seed() constructs a state-string using the "minimal" state-string format.
			 *  A state-string is a human-readable, ASCII representation of the 
			 *  PRNG's internal state. The \em minimal state-string mimics 
			 *  GNU's implementation of std::mt19937;
			 *  a single line containing space-separated words of state (unsigned integers),
			 *  with the number of words equalling \ref state_size,
			 *  and terminating with the \ref state_size itself (as a sanity check).
			 *  \code
					s_0  s_1  s_2 ... s_(state_size - 1) state_size  // Generic minimal state-string
				 \endcode
			 *  An example PRNG, using \c word_size = 64 and \c state_size = 3,
			 *  could use the following minimal state-string:
			 *  \code  
				 8602057372317241997 13337802638347508439 10520579822156476364 3
				 \endcode
			 * 	 
			 *  If the PRNG's state has more information, it's \c friend \c operator>> 
			 *  will have have to choose repeatable default behavior.
			 *  For example, Seed() does not supply xorshift1024_star its
			 *  state variable \a p, which activates the default 
			 *  \a p=0 choice of xorshift1024_star::operator>>.
			*/
			void Seed();

			/*! @brief Seed the PRNG from a state-string stored as the first line of an
			 *  ASCII or UTF-8 file (e.g. from WriteState()).
			 * 
			 *  \note A user-created, manual seed should mimic the minimal state-string format of Seed().
			 * 
			 *  \warning If prng_t has an full state-string which is longer than the minimal format (e.g. xorshift2014_star)
			 *  a user-created, manual seed \em should \em not try to supply this extra information.
			 *  The user should only supply the minimal state-string. 
			 *  It is the responsibility of the PRNG to repeatably configure itself 
			 *  given only the minimal state-string.
			 * 
			 *  @param filePath 
			 *  The full path to the ASCII file containing the state-string (on the first line).
			*/
			void Seed_FromFile(std::string const& filePath);
			
			/*! @brief Seed the PRNG from a state-string 
			 *  (e.g. from GetState()).
			 * 
			  *  \note A user-created, manual seed should mimic the minimal state-string format of Seed().
			 * 
			 *  \warning If prng_t has an full state-string which is longer than the minimal format (e.g. xorshift2014_star)
			 *  a user-created, manual seed \em should \em not try to supply this extra information.
			 *  The user should only supply the minimal state-string. 
			 *  It is the responsibility of the PRNG to repeatably configure itself 
			 *  given only the minimal state-string.
			 * 
			 *  @param stateString
			 *  The state-string.
			*/
			void Seed_FromString(std::string const& stateString);
						
			/*! @brief Write the state-string of the PRNG to the first line of an ASCII file
			 *  (e.g. for future reseeding by Seed_FromFile()).
			 * 
			 *  \warning Overwrites the file without warning. 
			 *  Will not attempt to create missing directories.
			 * 
			 *  \throws Throws std::ifstream::failure if the file cannot be opened
			 * 
			 *  @param filePath 
			 *  The full path to the file.
			*/
			void WriteState(std::string const& filePath);
			
			/*! @brief Return the full ASCII state-string of the PRNG
			 *  (e.g. for future reseeding by Seed_FromString()).
			*/
			std::string GetState();
			
			#if PRNG_CAN_JUMP
			/*! @brief Return a vector of PRNG state-strings (see GetState()) 
			 *  for use by a cadre of threads running independent generators.
			 *  
			 *  To give each thread its own unique generator, 
			 *  we can repeatedly Jump() a seeded generator, then record its state.
			 *  This is a repeatable way to obtain a vector of orthogonal PRNG states. 
			 *  The states are returned as state-strings because this is
			 *  the easiest way to portably communicate the PRNG state to each thread
			 *  ( e.g. if it is launched on the cloud).
			 *  The first state-string in the vector is the original state of the PRNG,
			 *  and each additional state-string is Jump()-ed one additional time. 
			 *  The PRNG ends in a state not represented in the vector, 
			 *  so it can continue to be used without risk of future collision.
			 *  
			 *  Storing the original state of the generator allows the 
			 *  same vector of state-strings to be regenerated in the future
			 *  (because Jump() is repeatable). Simple re-seed from the 
			 *  original state, then call this functions again.
			 *  
			 *  \param numThreads
			 *  The number of Jump()-ed state-strings to return.
			 * 
			 *  \return A vector of state-strings separated by one call to Jump().
			 *  The first state-string is the original state of the PRNG.
			*/
			std::vector<std::string> GetState_JumpVec(size_t const numThreads)
			{
				std::vector<std::string> stateVec;
				
				for(size_t i = 0; i < numThreads; ++i)
				{
					stateVec.push_back(GetState());
					this->Jump();
				}
				
				return stateVec;	
			}
			#endif
	};

	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	/*! @brief An implementation of the xorshift1024* 
	 *  64-bit pseudo-random number generator (PRNG).
	 * 
	 *  @author Sebastiano Vigna (vigna@acm.org)
	 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
	 *
	 *  KP's implementation of SV's 2017 edit of 
	 *  <a href="http://xoroshiro.di.unimi.it/xorshift1024star.c">
	 *  xorshift1024* </a> (period \f$ p = 2^{1024} - 1\f$),
	 *  with small changes to the code syntax and an API that mimics std::mt19937.
	 * 
	 *  \attention Per SV:
	 *  "... the two lowest bits of this generator are LFSRs of degree
	 *  1024, and thus will fail binary rank tests. The other bits needs a much
	 *  higher degree to be represented as LFSRs. We suggest to use a 
	 *  sign test to extract a random Boolean value, and
	 *  right shifts to extract subsets of bits..".
	 *  Hence, KP conservatively chooses not use the last 2 bits inside pqRand::engine.
	 * 
	 *  \attention The state must be seeded so that it is not everywhere zero, 
	 *  otherwise it will always return 0 (a property it shares with \c std::mt19937).
	 *  
	 *  SV claims that the generator has good quality and speed,
	 *  passing TestU01 with <a href="http://xoroshiro.di.unimi.it/#quality">
	 *  less errors than MT19937</a>.
	 *  KP independently tested this implementation using dieharder 3.31.1
	 *  \code
			 ./xorshift1024star_tester.x | dieharder -g 200 -a -k 2
		 \endcode
	 *  and it passed with flying colors (no failures).
	 *  D. Lemire reports 
	 *  <a href="https://lemire.me/blog/2017/09/15/the-xorshift1024-random-number-generator-fails-bigcrush/">
	 *  on his blog </a> that xorshift1024star fails Big Crush. 
	 *  SV implies <a href="http://xoroshiro.di.unimi.it/#FUD"> on his webpage </a> that 
	 *  this is the same issue he reported with the last two bits.
	 *  This potential failure is definitely something to watch.
	*/	 
	class xorshift1024_star
	{
		/* No doubt an SIMD version of this generator would be faster for massively parallel sampling.
		 * This would produce 16 new random words in one large batch, bu would require doubling the state
		 * 	s0 = current state
		 * 	s1 = last state
		 * We seed the generator as normal
		 * 	s0 = seed
		 * 	s1 = s0 (rotated left one index) ... for later use
		 * To warm up, we run the generator 16 times (using state = s0 and the original algorithm).
		 * This makes sure s0 has incremented in every index. Now every time we run the generator,
		 * we replace s0 -> s0[i] and s1 -> s1[i], then loop i = [0,16), which should be auto-vectorized.
		 * Afterwards, s0 is unchanged, s1 stores the new state, which we copy to s2 for return.
		 * Now s1 = s0 (rotated left one index) and s0 = s2 (less copying by swapping pointers, but maybe not faster).
		 * I'm pretty sure this is what Lemire did in (https://github.com/lemire/SIMDxorshift),
		 * and he reported a doubling of the sampling speed. I would expect more than doubling, 
		 * given the 4 times larger state and increased vectorized advantage.
		 * However, speeding up the PRNG will just further highlight that
		 * the real bottleneck is transforming the uniform integers into useful numbers,
		 * because it takes a lot of cycles and generally involves
		 * branching logic in the various math functions (i.e. I'm not aware of
		 * good, universal SIMD exp, log, pow with 1 ULP or error).
		*/ 		
		public:
			typedef uint64_t result_type; //!< @brief The unsigned integer type returned by the generator
			size_t static constexpr word_size = 64; //!< @brief Number of \em bits per PRNG \em word (i.e. result_type)
			size_t static constexpr state_size = 16; //!< @brief Number of \em words in the PRNG state (rounded up).
			
		 private:
			// Lots of deep magic in this class, not much to comment for SV's code			
			
			// Use std::array (versus C array) to give xorshift1024_star
			// an implicitly defined assignment (operator =),
			// while remaining compatible with C-style array syntax in SV's code
			std::array<uint64_t, state_size> state;
			uint64_t p;
			
		public:
			/*! @brief Power up the generator, but leave the state un-initialized.
			 *  
			 *  This usually places the generator in valid, but undefined, state.
			 *  The user \em must subsequently seed the generator using operator>>().
			 *  
			 *  \warning Don't forget to seed the generator.
			*/
			xorshift1024_star():p(0) {}
			// We must have p < state_size, so make sure to initialize it.
			// That way, an unintentional call to an un-seeded generator won't cause a seg fault.
			
			virtual ~xorshift1024_star() {}			
			
			//! @brief The smallest value this PRNG can return
			static constexpr size_t min() {return std::numeric_limits<result_type>::min();}			
			//! @brief The largest value this PRNG can return		
			static constexpr size_t max() {return std::numeric_limits<result_type>::max();}
			
			uint64_t operator()(); //!< @brief Return the next 64-bit, unsigned integer
			
			/*! @brief Quickly jump the state of the generator forward by \f$ 2^{512} \f$ calls.
			 * 
			 *  This allows one to generate \f$ 2^{512} \f$ parallel instances of the generator
			 *  with sequences that won't collide until they are called \f$ 2^{512} \f$ times.
			 *  This is useful when many instances of xorshift1024_star are running in 
			 *  parallel threads (see \ref jumping_PRNG).
			 *   
			 *  \note Obviously, the state of the generator is jumped 
			 *  without actually calling it \f$ 2^{512} \f$ times 
			 *  (since, at a billion calls per second, that would take
			 *  \f$ 4\times10^{137} \f$ years ... a bit too long).
			 * 
			 *  \attention This method seems like much deeper magic than the generator itself,
			 *  so KP cannot guarantee that it actually works as advertised by SV.
			 *  However, after the <a href="http://stackoverflow.com/questions/34574701/xorshift1024-jump-not-commutative">
			 *  Jan 2016 patch </a>,
			 *  it is definitely commutative 
			 *  [(call, Jump, call) == (Jump, call, call)],
			 *  which certainly implies that it works.
			*/ 
			void Jump();
			
			// Declare stream operators as friends, so they can access 
			// the private members. This only declares that they are friends, 
			// but does not actually declare the functions themselves.		
			friend std::ostream& operator << (std::ostream& stream, xorshift1024_star const& gen);
			friend std::istream& operator >> (std::istream& stream, xorshift1024_star& gen);	
	};
		
	// Actually declare the friend functions for xorshift1024_star
	// (previously declared as friends, but not as actual functions).
	
	/*! @brief Output the complete human-readable state-string of the generator to a stream.
	 * 
	 *  The full state of xorshift1024_star adds the index iterator \a p
	 *  to the the minimal state-string (as described in seeded_uPRNG::Seed):
	 *  \code
			s_0  s_1  s_2 ... s_15  16  p                     // Complete state-string for xorshift1024*
		 \endcode
	 *  \a p is determined from the number of times the generator has been called (modulo 16),
	 *  not from a random source of entropy, which is why it doesn't count towards \ref state_size.
	 * 
	 *  @param stream 
	 *  The stream to which the state-string is written.
	 *  
	 *  @param gen
	 *  The generator whose state is written.
	*/
	std::ostream& operator << (std::ostream& stream, xorshift1024_star const& gen);
	
	/*! @brief Seed the state of the generator from a stream, 
	 *  assuming the state-string format of operator<<().
	 *  
	 *  If \a p is missing from the state-string (e.g. if it is the minimal state-string),
	 *  a default value of \a p = 0 is used.
	 *  
	 *  \throws Throws seed_error if the stream runs out before the state is full,
	 *  or if the state-string is otherwise malformed.
	 * 
	 *  \warning Even though the state should have high entropy, 
	 *  and should definitely not be all zeroes (or all ones),
	 *  this function \em only verifies that a supplied state-string 
	 *  can be parsed, and never that it is good or even valid.
	 *  So don't supply a bad seed.
	 * 
	 *  @param stream 
	 *  The stream from which the state-string is read.
	 *  
	 *  @param gen
	 *  The generator whose state is seeded.
	*/
	std::istream& operator >> (std::istream& stream, xorshift1024_star& gen);	
		
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
			
	/*! @brief The random number engine used by all distributions in pqRand.
	 * 
	 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
	 * 
	 *  \ref engine wraps the seeded_uPRNG in an API designed for the 
	 *  improved inversion method. This provides access to both 
	 *  uneven and even samples from \f$ U(0,1) \f$, 
	 *  a random bool that doesn't waste bits, and a simple way to apply a random sign.
	 *  For the most part, these tools are used internally by pqRand's distributions, 
	 *  but can also be used to construct custom sampling objects.
	 * 
	 *  Arbitrary states of pqRand_engine can be stored using WriteState() 
	 *  and re-seeded using Seed_FromFile(). See seeded_uPRNG for more details.
	 * 
	 *  \warning Future changes to the API are not anticipated, but may still be possible
	*/ 
	class engine : public seeded_uPRNG<PRNG_t>
	{
		public:
			//! @brief The lowest few bits may have linear dependencies.
			result_type static constexpr badBits = 2;
		
		private:
			result_type static constexpr numBitsPRNG = word_size;
			result_type static constexpr numBitsMantissa = std::numeric_limits<real_t>::digits;
			
			static_assert((numBitsPRNG >= numBitsMantissa), 
				"pqRand::PRNG_t must supply as many bits as the mantissa of pqRand::rand_t can hold.");
				
			static_assert(((PRNG_t::min() == std::numeric_limits<result_type>::min())
				and ((PRNG_t::max() == std::numeric_limits<result_type>::max())
					or ((std::numeric_limits<result_type>::digits > PRNG_t::word_size)
					and (PRNG_t::max() == (result_type(1) << PRNG_t::word_size) - 1)))),
				"pqRand::PRNG_t must fill all of its digits.");
			
			result_type static constexpr bitShiftRight_even = 
				numBitsPRNG - numBitsMantissa;
					
			// For the next 2, must convert to signed int to make negative,
			// because -uint => huge uint
			real_t static constexpr scaleToU_even = 
				real_t(std::exp2(-int(numBitsMantissa))); 
												
			real_t static constexpr scaleToU_uneven = 
				real_t(std::exp2(-int(numBitsPRNG)));
				
			result_type static constexpr replenishBitCache = (result_type(1) << (badBits - 1));
			// NOTE: if badBits == 0, replenishBitCache == 0 (a shift left of -1).
			// This may output a warning, but should still be valid.
						
			// We need to fill the mantissa, with 2 bits for rounding:
			// a buffer bit and a sticky bit. The buffer bit must use 
			// a good bit of entropy, but the sticky bit is always set to 1, 
			// so we can use the last bad bit as the sticky bit
			result_type static constexpr numBitsOfEntropyRequired = 
				numBitsMantissa + 1 + ((badBits > 0) ? (badBits - 1): 1);
			// When the random number is less than minEntropy, we need more entropy.
			result_type static constexpr minEntropy = 
				(result_type(1) << (numBitsOfEntropyRequired - 1));
						
			result_type bitCache; //! A cache of random bits for RandBool
			result_type cacheMask; // Selects one bit from bitCache, for RandBool
		
			// Top up the entropy when randUint does not have enough for an uneven variate
			real_t U_uneven_TopUpEntropy(result_type randUint);
			
			// Redefine the base class virtuals, because we need to 
			// store/refresh the state of the bitCache when we write/seed
			virtual void Seed_FromStream(std::istream& stream);
			virtual void WriteState_ToStream(std::ostream& stream);
			
			// We must always default-initialize the bitCache in the same way
			void DefaultInitializeBitCache();
			
		public:
			/*! @brief Construct the engine; auto-seed if requested.
			 * 
			 *  To auto-seed \ref engine, construct like:
			 *  \code
					pqRand::engine gen1();
					pqRand::engine gen2(true);
				 \endcode
				 To supply a seed, construct like:
			 *  \code
					pqRand::engine gen3(false); // Don't auto-seed (leaving gen3 is an undefined state).
					gen3.Seed_FromFile("./archive/Jan2017/run123.seed"); // Seed from a file (gen3 is no longer undefined).
				 \endcode
			 * 
			 *  \note If one intends to supply a user-generated seed
			 *  (as opposed to storing/reusing an auto-seed),
			 *  they should use the \em minimal state-string format
			 *  (which is described in Seed()).
			 *  
			 *  @warning \ref engine contains more state information than the 
			 *  minimal state-string. This extra information should only 
			 *  be output by WriteState(), and never generated by the user.
			 * 
			 *  @param autoSeed 
			 *  Perform an autoSeed (\c true) by calling Seed()
			 *  or defer seeding till later (\c false).
			*/
			engine(bool const autoSeed = true):
				seeded_uPRNG(autoSeed)
				// Let the base class do the autoSeed, because although
				// engine redefines Seed_FromStream, which the base class ctor cannot access,
				// Seed is only passing the minimal state-string, 
			{
				// We must handle everything not seeded by the super-class,
				// which seeded the PRNG from the minimal state-string
				DefaultInitializeBitCache();				
			}
			
			virtual ~engine() {}
			
			/*! @brief Return the result of an ideal coin flip
			 * 
			 *  This uses the PRNG efficiently, using 1 bit of randomness per \c bool.
			*/
			bool RandBool(); 
			
			/*! @brief Assign the victim a random sign (+/-) (by \em reference) using RandBool().
			 * 
			 * \param victim 
			 * The floating-point number which is assigned a random sign (+/-).
			 * 
			 * \return the altered victim
			*/ 			
			real_t ApplyRandomSign(real_t& victim)
			{
				if(RandBool()) victim = -victim;
				return victim;
			}
			
			/*! @brief Assign the victim a random sign (+/-) using RandBool().
			 * 
			 * \param victim 
			 * The floating-point number which is assigned a random sign (+/-).
			 * 
			 * \return the altered victim
			*/ 
			real_t ApplyRandomSign(real_t&& victim)
			{
				ApplyRandomSign(victim);
				return victim;
			}
			
			// The uniform variate functions are MUCH faster if they are defined in the header, 
			// because the compiler can do more optimization depending on their use.
			// The TopUpEntropy function is rare, so leave it to the library to reduce inline-ing
							
			/*! @brief Draw an uneven uniform variate from \f$ U(0, 1] \f$.
			 * 
			 *  Draw a random \em real number from \f$ U(0, 1] \f$, 
			 *  then round to the nearest floating point.
			 *  
			 *  \note 1 is half as probable as its next-door neighbor
			*/
			real_t U_uneven() 
			{
				result_type const randUint = (*this)();

				if(randUint < minEntropy) // If randUint lacks enough entropy, top it up
					return U_uneven_TopUpEntropy(randUint);
				else 
					return scaleToU_uneven * real_t(randUint bitor result_type(1));
			}
			
			/*! @brief Draw an uneven uniform variate from \f$ U(0, 0.5] \f$.
			 * 
			 *  Draw a random \em real number from \f$ U(0, 0.5] \f$, 
			 *  then round to the nearest floating point.
			 *  
			 *  \note 0.5 is half as probable as its next-door neighbor
			*/
			real_t HalfU_uneven() 
			{
				result_type const randUint = (*this)();

				if(randUint < minEntropy) // If randUint lacks enough entropy, top it up
					return real_t(0.5) * U_uneven_TopUpEntropy(randUint);
				else 
					return real_t(0.5) * scaleToU_uneven * real_t(randUint bitor result_type(1));
			}
			
			/*! @brief Draw an even uniform variate from \f$ U[0, 1)\f$
			 * 
			 *  Partition \f$ U[0, 1) \f$ in increments of machine \f$ \epsilon \f$
			 *  and draw uniformly from this sample-space.
			*/
			real_t U_even()
			{
				return scaleToU_even * real_t((*this)() >> bitShiftRight_even);
			}
	};	
};

#endif
