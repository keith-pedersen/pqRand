#include "pqRand.hpp"
#include <stdexcept>
#include <sstream>
#include <iostream>

using namespace pqr;

int main()
{	
	// This package uses C++11.
	
	// The pqRand_engine is the main module needed by every distribution.
	// It does the random integer generation, converts integers to U_Q, 
	// and does the random coin flips needed by the quantile flip-flop.
	// All pqr::distributions require a pqRand_engine to be supplied by 
	// reference when drawing from them. This is a nearly identical API to 
	// the std::distributions of C++11.
	
	// The pqRand_engine forces an initial seed AS LARGE as the generators state
	// (e.g. a 32-bit seed doesn't fill up 1024 bits of state).
	// If no seed is supplied, it uses std::random_device, 
	// (which uses /dev/urandom on Linux) to get actual random entropy.
	// This is a standardized way to get a good source of initial entropy;
	// if a source of true entropy is unavailable, the implementation is 
	// supposed to implement a decent PRNG for you.
	// A seed from a file or an std::istream can also be supplied, 
	// provided they are in the right format (see pqRand.hpp).
	// This allows the user to supply a custom seed.
	
	// Let's look at a few ways to start up a pqRand_engine
	
	{	
		// No args ... seed gen1 using std::random_device,
		pqRand_engine gen1;
		
		// Store gen1's seeded initial state to a file, for re-use;
		gen1.WriteState("test.seed");
		
		// Seed another generator from the stored seed
		pqRand_engine gen2("test.seed");
		
		// We can also store a seed in any std::ostream. 
		// A stringstream is a less permanent place to store a seed
		std::stringstream seed;
		gen2.WriteState(seed);

		// And we can seed yet another generator from gen2;
		pqRand_engine gen3(seed);
		
		// We can also directly copy generators
		pqRand_engine gen4 = gen3;
		
		printf("\n Seed test\n");
		printf("--------------------------------------------------------------------------------\n");
		printf("          gen1                  gen2                  gen3                  gen4\n");
		for(size_t i = 0; i < 5; ++i)
		{
			printf("%lu  %20lu  %20lu  %20lu  %20lu\n", i, gen1(), gen2(), gen3(), gen4());
		}
		printf("\n");
	}
	
	// We may also wish to have a bunch of parallel threads.
	// But if every thread calls the same generator, the determinism of 
	// the PRNG is broken, because thread scheduling is not deterministic.
	// The solution is to give each thread its own generator.
	// the pqRand package uses xorshift1024* (see pqRand.hpp)
	// because its state size (1 kB) is much smaller than common generators
	// like the Mersenne Twisters (20 kB). 
	// However, choosing random seeds for each thread might cause 
	// accidental collisions in their sequences.
	// This is why the author of xorshift1024* created the jump function, 
	// which advances the generator by 2^512 calls 
	// (obviously without calling that many times).
	// Here is an example of creating 8 generators,
	// each separated by a single jump.
	
	{
		std::vector<pqRand_engine> engineList;
		
		engineList.emplace_back("test.seed");
		
		for(size_t i = 0; i < 4; ++i)
		{
			engineList.emplace_back(engineList.back());
			engineList.back().Jump();
		}
		
		printf("\n Jump test\n");
		printf("--------------------------------------------------------------------------------\n");
		printf("    3 calls, then jump all generators to the same state, so their back in sync\n\n");
		
		for(size_t i = 0; i < 3; ++i)
		{
			printf("%10lu", i);
			for(auto& gen : engineList)
			{
				printf("	 %20lu", gen());
			}
			printf("\n");
		}
		
		printf("jump-sync ");
		for(size_t i = 0; i < engineList.size(); ++i)
		{
			engineList[i].Jump(engineList.size() - i - 1);
			printf("	 %20lu", engineList[i]());			
		}
		printf("\n\n");
	}
	
	printf("\n\n Utilities\n");
	printf("--------------------------------------------------------------------------------\n");
	printf("    pqRand_engine gives access to uint64_t, U_Q, HalfU_Q, and random bool (as well as U_S)\n");
	
	pqRand_engine gen;
	// The pqRand_engine::operator()() gives us direct access to 
	// the uint64_t calculated by mt19937_64
	// To get (0., 1.) from pqRand_engine, we use the function U_Q(). 
	// HalfU_Q() gives a number in (0, 0.5] (for use in certain flip-flops).
	printf("gen():          %lu\n", gen());
	printf("gen.U_Q():      %.17e\n", gen.U_Q());
	printf("gen.HalfU_Q():  %.17e\n", gen.HalfU_Q());
	printf("gen.RandBool(): %d  %d  %d  %d  %d  %d\n", 
		gen.RandBool(), gen.RandBool(), gen.RandBool(), gen.RandBool(), gen.RandBool(), gen.RandBool());
	// There is also a function to apply a random sign to a double passed by refernce (ApplyRandomSign())
	
	printf("\n\n Distributions\n");
	printf("--------------------------------------------------------------------------------\n");
	printf("    pqRand has built-in access to a number of distribution\n");
			
	// Now let's play with some distributions.
	// In all cases, I used the explicit PDFs given by Wikipedia, 
	// which may lead to some issues of convention.
	standard_normal stand; 				// mu = 0, sigma = 1, always
	normal norm(-1.5, 3.1); 			// mu = -1.5, sigma = 3.1
	log_normal logNorm(2.71, 0.66);  // mu = 2.71, sigma = 0.66
	weibull weib(4.56, 1.23); 			// lambda = 4.56, k = 1.23
	pareto par(3.33, 4.); 				// x_m = 3.33, alpha = 10.
	
	// I have a histogram class which I use to bin these, 
	// but I figure you probably have your own preferred way of looking at them
	printf("\n");
	printf(" std_norm      normal      logNormal     weibull      pareto\n");
	for(int i = 0; i < 10; ++i)
		printf("% .3e   % .3e   % .3e   % .3e   % .3e\n", 
		stand(gen), norm(gen), logNorm(gen), weib(gen), par(gen));
	printf("\n\n");
}
