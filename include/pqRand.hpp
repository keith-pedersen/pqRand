/*!
 *  @file pqRand.hpp
 *  @brief Establishes pqRand and defines the PRNG engine needed by distributions
 * 
 *  The important methods implemented here are the "quasiuniform sample" of 
 *  the semi-inclusive unit interval (pqRand::engine::U_Q) and
 *  the random coin-flip (pqRand::engine::RandBool).
 * 
 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
 *  @date 2017
*/
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
 *  @version 0.4.3 (beta); June 2017
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
 *  This package is designed to improve sampling the following distributions:
 *     - \ref pqRand::normal "normal"
 *     - \ref pqRand::log_normal "log-normal"
 *     - \ref pqRand::exponential "exponential"
 *     - \ref pqRand::weibull "weibull"
 *     - \ref pqRand::pareto "pareto"
 *    
 *  @anchor theory
 *  The principle theory behind pqRand is presented in the paper
 *  "Conditioning your quantile function," https://arxiv.org/abs/1704.07949
 *  (distributed in doc/).
 *  pqRand permits a discrete, high-precision sample from a 
 *  selection of distributions which possess an analytic quantile function 
 *  (or, as is the case for the normal and log-normal, 
 *  which can be drawn using a quantile function of an underlying distribution). 
 *  pqRand uses the "quantile flip-flop", as described in the paper. 
 *  For quality control purposes, the user is "forced" to use the 
 *  built-in pseudo-random number generator (PRNG),
 *  although this can be changed (within limits, see pqRand::PRNG_t).
 * 
 * 
 *  Usage
 *  --------------------------------------------------------------------
 *  1. Create a pqRand::engine (one per thread) to draw random numbers.
 *     - "Auto-seed" is highly recommended (see pqRand::engine::Seed).
 *     - For auditing, only the first thread's seed needs to be stored (see pqRand::engine::Jump).
 *  2. Create as many distributions as you want.
 *  3. Call the distributions like functions, passing the engine as the argument.
 
    \code
		gen = pqRand::engine(); // Create the random engine with auto-seed (default constructor)
		dist = pqRand::pareto(2); // Instantiate the distribution
		
		std::vector<double> sample;
		for(size_t i = 0; i < 1024; ++i)
			sample.push_back(dist(gen)); // Repeatedly sample from the distribution
		
    \endcode
 *  For more examples, see "examples/".
 * 
 * 
 *  Python
 *  -------------------------------------------------------------------
 *  The Python/Cython API is 99% identical to the C++ API, so I have not 
 *  built a parallel documentation for Python. The two minor differences are 
 *  due to limitations of the Python language, and are enumerated at the top of "pYqRand.pyx".
 * 
 * 
 *  Features to consider in future versions
 *  --------------------------------------------------------------------
 *  - If a regular sample of \f$ U \f$ produces numbers with \f$ weight = 1 \f$, 
 *    supply a sample of \f$ U \f$ from the tail \f$(u \to 0)\f$,
 *    given a supplied \f$ weight \ll 1 \f$.
 * 
 * 
 *  Known bugs
 *  --------------------------------------------------------------------
 *  - None at this time.
*/
 
#ifndef PQ_RAND
#define PQ_RAND

#include <limits> // numeric_limits
#include <string>
#include <array>
#include <random> // random_device, mt19937

namespace pqRand //! @brief The namespace of the pqRand package
{
	/*! @brief The PRNG used by pqRand::engine
	 * 
	 *  The authors of pqRand have carefully chosen xorshift1024_star,
	 *  a PRNG with a sufficiently large state (but not too large), 
	 *  excellent uniformity properties, rigorously tested by others and 
	 *  vetted by the authors using the dieharder battery of tests.
	 *  
	 *  The motivation of the pqRand package essentially requires that
	 *  simple freedom to chose a PRNG be removed from the user, 
	 *  since such freedom would permit them to unwittingly destroy the 
	 *  benefits of pqRand. Nonetheless, it is still possible to 
	 *  override this decision by redefining PRNG_t, IFF the substituted class 
	 *  uses satisfies the \ref prng_requirements "\c prng_t requirements" 
	 *  of seeded_uPRNG.
	 * 
	 *  As such, \c std::mt19937_64 (the model for the xorshift102_star API)
	 *  will work as a drop-in replacement. But while 
	 *  \c std::mt19937_64 has good statistics, 
	 *  its state is too large for massive multi-parallelism.
	 *  Nonetheless, it is fully compatible with pqRand
	 *  (though not necessarily with pqRand_Example.cpp, which uses Jump()).
	 *  To switch to std::mt19937_64, simply change this typedef
	*/ 
	// typedef xorshift1024_star PRNG_t;
	// typedef std::mt19937_64 PRNG_t;
	typedef std::mt19937 PRNG_t; // valid if real_t == float
	
	//! @brief Defined to allow internal testing with float/binary32.
	// typedef double real_t;
	typedef float real_t;
	
	/*! @brief A wrapper for a PRNG (32 or 64 bit), providing a seeding interface (used by \ref engine).
	 * 
	 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
	 * 
	 *  The seeding functions fill the \em entire state of the PRNG
	 *  (because you can't properly seed a PRNG with 128 bytes of state using a 4-byte int).
	 *  The recommended method is using Seed() to "auto-seed" the PRNG.
	 * 
	 *  The state of the PRNG can be converted to an ASCII "state-string" by GetState().
	 *  This state-string can be written directly to a file by WriteState(). 
	 *  This permits reseeding by Seed_FromString() or Seed_FromFile() (respectively). 
	 *  The PRNG can also be seeded from a user-generated state-string, 
	 *  provided that it is in the correct format (see Seed()).
	 * 
	 *  \warning Seed_FromFile() and Seed_FromString() do not check the format of
	 *  state-strings before forwarding them to \c prng_t.
	 *  
	 *  \anchor SeedWrite_Virtual
	 *  Seed_FromStream() and WriteState_ToStream() are
	 *  simple wrappers to \c prng_t 's \c operator<< and \c operator>>.
	 *  They are the actual worker functions called by the public Seed/State functions, 
	 *  and are \b virtual so that a derived PRNG class may add 
	 *  more information to its internal state (beyond the state of the underlying PRNG)
	 *  without redefining the public Seed/State functions.
	 * 
	 *  \internal There are some extra tabs in here to get the list-formatting right in Doxygen. 
	 *            They seem to be needed for the bullets that contain the code examples. \endinternal
	 *  @anchor prng_requirements
	 *  @param prng_t    
	 *  A PRNG class supplying \c result_type. It must:
	 *  	- have a default constructor 
	 *    	- possess the following static fields
	 * 	\code
	       static const word_size; // Number of bits per PRNG word
	       static const state_size; // Number of words of seed entropy required.
	       typedef result_type; // The unsigned integer type returned by the generator
			\endcode
	 *    	- implement the following functions to 
	 *      input/output the state of the generator in ASCII format
	 *  \code
	       friend operator>>(std::istream&, prng_t&);
	       friend operator<<(std::ostream&, prng_t&);	       
	    \endcode
	 *		- support the "minimal" format used by Seed().
	*/
	template<class prng_t>
	class seeded_uPRNG : public prng_t
	{
		public: 
			using prng_t::word_size;
			using prng_t::state_size;
			using typename prng_t::result_type;		
		
			static_assert(((word_size == 32) or (word_size == 64)), "pqRand::seeded_uPRNG requires either a 32-bit or 64-bit prng_t");
			static_assert((std::numeric_limits<result_type>::digits >= word_size), "pqRand::seeded_uPRNG: prng_t::result_type must be able to accomadate a full word");
		
		protected:
			/*! @brief Set the state of the PRNG from a state-string stored in a stream, 
			 *  using \c prng_t's \c operator>>.
			 * 
			 *  See \ref SeedWrite_Virtual "the detailed class description".
			 * 
			 *  @param stream 
			 *  The state input.
			*/
			virtual void Seed_FromStream(std::istream& stream);
			
			/*! @brief Write the state-string of the PRNG to \a stream, using \c prng_t's \c operator<<.
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
			 *  	(see engine::engine for an example).
			 * 
			 *  @param autoSeed 
			 *  Perform an autoSeed (\c true) by calling Seed()
			 *  or defer seeding till later (\c false).
			*/
			explicit seeded_uPRNG(bool const autoSeed = true):
				prng_t()
			{
				if(autoSeed) Seed(); // else defer seed to user
			}
			// NOTE: Access to non-auto Seed from inside the ctor is difficult.
			// They use the same arguments, so which do you want?
			// Plus, Cython can't overload functions, so we need a uniform ctor.
			
			/*! @brief Auto-seed the generator using \c std::random_device.
			 * 
			 *   <a href="http://en.cppreference.com/w/cpp/numeric/random/random_device">
			 *  \c std::random_device </a> 
			 *  is intended to provide access to a small amount of \em true randomness
			 *  (or at least a non-deterministic PRNG).
			 *  On many GNU/Linux systems, \c std::random_device uses \c \\dev\\urandom.
			 *  Seed() calls \c std::random_device until 
			 *  the generator's state is filled.
			 * 
			 *  Seed() constructs its state-string using the "minimal" format 
			 *  described in xorshift1024_star::operator<<. Thus, 
			 *  for \ref prng_t = xorshift1024_star, 
			 *  the state variable \a p is not supplied,
			 *  activating the default \a p behavior of xorshift1024_star::operator>>.
			*/ 
			void Seed();
	
			/*! @brief Seed the PRNG from a state-string stored in an ASCII file 
			 *  (e.g. from WriteState()).
			 * 
			 *  @note To create a user-defined seed, mimic the state-string format of Seed().
			 * 
			 *  @param fileName 
			 *  The full path to the ASCII file containing the state-string (on the first line).
			*/
			void Seed_FromFile(std::string const& fileName);
			
			/*! @brief Seed the PRNG from a state-string 
			 *  (e.g. from GetState()).
			 * 
			 *  @note To create a user-defined seed, mimic the state-string format of Seed().
			 * 
			 *  @param seedString
			 *  The state-string.
			*/
			void Seed_FromString(std::string const& seedString);
						
			/*! @brief Write the state-string of the PRNG to an ASCII file
			 *  (e.g. for future reseeding by Seed_FromFile()).
			 * 
			 *  \warning Overwrites the file without warning. 
			 *  Will not attempt to create missing directories.
			 *  Will throw an exception if the file cannot be opened.
			 * 
			 *  @param fileName 
			 *  The full path to the file.
			*/
			void WriteState(std::string const& fileName);
			
			/*! @brief Return the ASCII state-string of the PRNG
			 *  (e.g. for future reseeding by Seed_FromString()).
			*/
			std::string GetState();
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief An implementation of the xorshift1024* 
	 *  64-bit pseudo-random number generator (PRNG).
	 * 
	 *  @author Sebastiano Vigna (vigna@acm.org)
	 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
	 *
	 *  KP's implementation of SV's <a href="http://xoroshiro.di.unimi.it/xorshift1024star.c">
	 *  xorshift1024* </a> (period \f$ p = 2^{1024} - 1\f$),
	 *  with small changes to the code syntax and a new API that mimics std::mt19937.
	 *  
	 *  SV claims that the generator has good quality and speed,
	 *  passing TestU01 with <a href="http://xoroshiro.di.unimi.it/#quality">
	 *  less errors than MT19937</a>.
	 *  KP independently tested this implementation using dieharder 3.31.1
	 *  \code
	       ./xorshift2014star_tester.x | dieharder -g 200 -a -k 2
	    \endcode
	 *  and it passed with flying colors (no failures).
	 * 
	 *  \attention Per SV:
	 *  "... The three lowest bits of this generator are LSFRs, and thus
	 *  they are slightly less random than the other bits. 
	 *  We suggest to use a sign test to extract a random Boolean value."
	 * 
	 *  The last 3 bits also passed dieharder
	 *  (concatenating the last-3 bits from three sequential calls)
	 *  \code
	       ./xorshift2014star_last3bits_tester.x | dieharder -g 200 -a
	    \endcode
	 *  so perhaps they're not as bad as SV implies.
	 *  But for now, we conservatively do not use the last 3 bits inside pqRand::engine.
	*/	 
	class xorshift1024_star
	{
		public:
			size_t static constexpr word_size = 64; //!< Number of \em bits per PRNG \em word
			size_t static constexpr state_size = 16; //!< Number of \em words of seed entropy required.
			typedef uint64_t result_type;
		
		 private:
			// Lots of deep magic in this class, not much to comment for SV's code			
			uint64_t static const JUMP[state_size];
			
			// Use std::array (versus C array) to utilize xorshift1024_star's
			// implicitly defined assignment (operator =),
			// while remaining compatible with C-style array syntax in SV's code
			std::array<uint64_t, state_size> state; 
			uint64_t p;
			
		public:
			/*! @brief Power up the generator, but leave the state un-initialized.
			 *  
			 *  This places the generator in a valid, but undefined, state.
			 *  The user must subsequently seed the generator using operator>>().
			 *  
			 *  \warning Don't forget to seed the generator.
			*/
			xorshift1024_star():p(0) {}
			// We must have p < state_size, so make sure to initialize it.
			// This is a bad state, but valid 
			// ... unless the uninitialized state happens to be all zeros.
						
			uint64_t operator()(); //!< @brief Return the next 64-bit, unsigned integer
			
			/*! @brief Quickly jump the state of the generator forward by \f$ 2^{512} \f$ calls.
			 * 
			 *  This is useful when many instances of xorshift1024_star are running in 
			 *  parallel threads, since feeding each its own random seed might accidentally cause
			 *  the sequences from one pair of generators to overlap, 
			 *  thus become correlated (combinatorially, this is 
			 *  too probable and too calamitous to ignore).
			 *  Alternatively, storing the seeds from \em every thread
			 *  in \em every simulation --- for auditing purposes --- 
			 *  would require a lot of storage.
			 * 
			 *  The Jump() method only requires storing one seed per simulation: 
			 *  \code
					size_t const numThreads = 16; // A simulation with 16 threads
					std::vector<pqRand::engine> engineVec; // A vector storing each threads generator
					
					engineVec.emplace_back(); // Auto-seed the first generator ...
					engineVec.back().WriteState("mySimulationID.seed"); // then store its seed.
					
					// Now create a bunch of orthogonal generators
					for(size_t threadIndex = 1; threadIndex < numThreads; ++threadIndex)
					{
						engineVec.emplace_back(engineVec.back()); // Each engine starts as a copy of the last engine ...
						engineVec.Jump(); // but is jumped one additional time
					}
			    \endcode 
			 *  
			 *  \note Obviously, the state of the generator is jumped 
			 *  without actually calling it \f$ 2^{512} \f$ times 
			 *  (since, at a billion calls per second, that would take
			 *  \f$ 4\times10^{137} \f$ years ... a bit too long).
			 * 
			 *  \attention This method seems like much deeper magic than the generator itself,
			 *  so KP cannot guarantee that it actually works as advertised by SV.
			 *  Howerver, after the  <a href="http://stackoverflow.com/questions/34574701/xorshift1024-jump-not-commutative">
			 *  Jan 2016 patch </a>,
			 *  it is definitely commutative 
			 *  [(call, Jump, call) == (Jump, call, call)],
			 *  which certainly implies that it works.
			 
			*/ 
			void Jump();
						
			void Jump(size_t nTimes); //!< @brief Call Jump() \a nTimes.
			
			// Declare stream operators as friends, so they can access 
			// the private members
			
			/*! @brief Output the complete state-string of the generator to a stream.
			 * 
			 *  The state-string is an ASCII representation of the 
			 *  PRNG's internal state. In pqRand, the \em minimal state-string mimics 
			 *  GNU's implementation of std::mt19937:
			 *  a single line containing space-separated words of state,
			 *  with the number of words equalling \ref state_size,
			 *  and terminating with the \ref state_size itself.
			 *  \code
					s_0  s_1  s_2 ... s_(state_size - 1)  state_size  // Generic minimal state-string
					s_0  s_1  s_2 ... s_15  16                        // Minimal state-string for xorshift1024*
				 \endcode
			 *  An example PRNG, using \c word_size = 64 and \c state_size = 3, 
			 *  could use the following minimal state-string:
			 *  \code  
				 8602057372317241997 13337802638347508439 10520579822156476364 3
				 \endcode
			 * 	 
			 *  A PRNG may contain additional components of its internal state 
			 *  which are not accounted for in \ref state_size.
			 *  For xorshift1024_star, the exact state includes an index iterator 
			 *  \a p (where \a p < \ref state_size). Thus, its complete state-string is:
			 *  \code
					s_0  s_1  s_2 ... s_15  16  p                     // Complete state-string for xorshift1024*
				 \endcode
			 *  Because \a p < \ref state_size, it should not be determined 
			 *  from a random source of entropy, but from the number of times
			 *  the generator is called (modulus 16). This is why it is 
			 *  not counted towards \ref state_size.
			 * 
			 *  @param stream 
			 *  The stream to which the state-string is written.
			 *  
			 *  @param gen
			 *  The generator whose state is written.
			*/
			friend std::ostream& operator << (std::ostream& stream, xorshift1024_star const& gen);
			
			/*! @brief Seed the state of the generator from a stream, 
			 *  assuming the state-string format of operator<<().
			 *  
			 *  If \a p is missing from the state-string, 
			 *  a default value of \a p = 0 is used (see seeded_uPRNG::Seed()).
			 *  
			 *  @warning Throws an exception if the stream runs out before the state is full,
			 *  or if the stream does not follow the expected format.
			 * 
			 *  \warning Per SV:
			 *  "The state must be seeded so that it is not everywhere zero."\n
			 *  Nonetheless, the \em only checks made by operator<<() on a supplied state-string 
			 *  is that it can be parsed in the expected format.
			 *  When generating a user-defined seed (instead of using an auto-seed),
			 *  take care not to supply a bad seed. 
			 * 
			 *  @param stream 
			 *  The stream from which the state-string is read.
			 *  
			 *  @param gen
			 *  The generator whose state is seeded.
			*/
			friend std::istream& operator >> (std::istream& stream, xorshift1024_star& gen);	
	};
		
	// Actually declare the friend functions for xorshift1024_star
	// (previously declared as friends, but not as actual functions).
	std::ostream& operator << (std::ostream& stream, xorshift1024_star const& gen);
	std::istream& operator >> (std::istream& stream, xorshift1024_star& gen);	
		
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
			
	/*! @brief The random number engine used by all distributions in pqRand.
	 * 
	 *  \ref engine wraps the seeded PRNG in an API designed for the quantile flip-flop.
	 *  This provides access to both quasuniform and superuniform samples from \f$ U \f$, 
	 *  a random bool that doesn't waste bits, and a simple way to apply a random sign.
	 *  For the most part, these tools are used internally by prRand's distributions.
	 * 
	 *  Arbitrary states of pqRand_engine can be stored using WriteState() 
	 *  and re-seeded using Seed_FromFile(). See seeded_uPRNG for more details.			 
	 * 
	 *  \note Future changes to the API are not anticipated, but may still be possible
	*/ 
	class engine : public seeded_uPRNG<PRNG_t>
	{
		public: 
			using seeded_uPRNG<PRNG_t>::word_size;
			using seeded_uPRNG<PRNG_t>::state_size;
			using seeded_uPRNG<PRNG_t>::result_type;
			
		protected:
			result_type static constexpr numBitsPRNG = word_size;
			result_type static constexpr numBitsMantissa = std::numeric_limits<real_t>::digits;
			
			static_assert((numBitsPRNG > numBitsMantissa), "prRand::PRNG_t must supply more bits than the mantissa of pqRand::rand_t can hold");
			
			result_type static constexpr bitShiftRight_Superuniform = numBitsPRNG - numBitsMantissa;
			result_type static constexpr maxMantissa_Superuniform = (result_type(1) << numBitsMantissa);
		
			// Scale the mantissa into the correct interval
			real_t static constexpr scaleToU_Quasiuniform = 
				real_t(std::exp2(-int(word_size))); // convert to int to allow negative
			
			real_t static constexpr scaleToHalfU_Quasiuniform = 
				real_t(0.5) * scaleToU_Quasiuniform;
				
			real_t static constexpr scaleToU_Superuniform = 
				// epsilon = 2**-(P-1), we want 2**-P
				real_t(0.5) * std::numeric_limits<real_t>::epsilon();
			
			real_t static constexpr scaleToHalfU_Superuniform = 
				real_t(0.5) * scaleToU_Superuniform;
					
			// The last three bits of xorshift1024_star are not as good to use, per S. Vigna
			result_type static constexpr badBits = 3;
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
		
			// Draw a random, quasiuniform mantissa from (0, 2**53],
			// with 2**53 being half as probable as 2**53-1
			real_t RandomMantissa_Quasiuniform();
			//~ // Draw a random, superuniform mantissa from (0, 2**53],
			//~ // with 2**53 being half as probable as 2**53-1
			//~ real_t RandomMantissa_Superuniform();
			// Draw a random, superuniform mantissa from (0, 2**53),
			// in the canonical manner (equi-probable sample space)
			real_t RandomMantissa_Superuniform();
			
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
			 *  If one intends to supply a user-generated seed
			 *  (as opposed to storing/reusing an auto-seed),
			 *  they should use the \em minimal state-string format
			 *  (which is described in Seed()) using a \ref state_size of 16.
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
				// Defer automatic seed from base class because we need to use 
				// virtual Seed_FromStream to properly seed, which the
				// base-class can't access from inside its ctor.
			{
				if(autoSeed)	
					Seed();
				else // Defer seeding till later, but place the bitCache in a valid state
					DefaultInitializeBitCache();
			}
			
			bool RandBool(); //!< @brief Return the result of an ideal coin flip
			
			/*! @brief Give the victim a random sign (+/-), using RandBool().
			 * 
			 * @param victim 
			 * The floating-point number which is assigned a random sign (+/-).
			*/ 			
			void ApplyRandomSign(real_t& victim); 
				
			/*! @brief Draw from quasiuniform \f$ U(0, 1] \f$.
			 * 
			 *  Draw a random \em real number from \f$ U(0, 1] \f$, then round to the 
			 *  nearest floating point.
			 *  
			 *  \note 1 is half as probable as its next-door neighbor
			*/
			inline real_t U_Q() 
			{
				return scaleToU_Quasiuniform * RandomMantissa_Quasiuniform();
			}
			
			/*! @brief Draw from quasiuniform \f$ U(0, 0.5] \f$.
			 * 
			 *  Draw a random \em real number from \f$ U(0, 0.5] \f$, then round to the 
			 *  nearest floating point.
			 *  
			 *  \note 0.5 is half as probable as its next-door neighbor
			*/
			inline real_t HalfU_Q() 
			{
				return scaleToHalfU_Quasiuniform * RandomMantissa_Quasiuniform();
			}
			
			//~ /*! @brief Draw from superuniform \f$ U(0, 1] \f$.
			 //~ * 
			 //~ *  Partition \f$ U(0, 1] \f$ in increments of machine \f$ \epsilon \f$
			 //~ *  and draw uniformly from this sample-space.
			 //~ *  
			 //~ *  \note 1 is half as probable as its next-door neighbor
			//~ */
			//~ inline real_t U_S() 
			//~ {
				//~ return scaleToU_Superuniform * RandomMantissa_Superuniform();
			//~ }
			
			//~ /*! @brief Draw from superuniform \f$ U(0, 0.5] \f$.
			 //~ * 
			 //~ *  Partition \f$ U(0, 0.5] \f$ in increments of \f$ \epsilon/2 \f$
			 //~ *  and draw uniformly from this sample-space.
			 //~ *  
			 //~ *  \note 0.5 is half as probable as its next-door neighbor
			//~ */
			//~ inline real_t HalfU_S()
			//~ {
				//~ return scaleToHalfU_Superuniform * RandomMantissa_Superuniform();
			//~ }
			
			/*! @brief Draw from superuniform \f$ U(0, 1)\f$ (in the canonical manner).
			 * 
			 *  Partition \f$ U(0, 1) \f$ in increments of machine \f$ \epsilon \f$
			 *  and draw uniformly from this sample-space.
			*/
			inline real_t U_S() 
			{
				return scaleToU_Superuniform * RandomMantissa_Superuniform();
			}
	};	
};

#endif
