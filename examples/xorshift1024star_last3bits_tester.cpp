#include "pqRand.hpp"
#include <fstream>
#include <iostream>

using namespace pqr;

int main()
{
	uPRNG_64_Seeder<xorshift1024_star> gen;
	gen.Seed();
	
	// Continously output random numbers to std::out, for use like ... 
	// 	./xorshift1024star_last3bits_tester.x | dieharder -g 200 -a
	// Preliminary results (i.e. running the above line of code)
	// show that the last 3 bits of xorshift1024* are not that bad	
	while(true)
	{
		// Make a mutant 9 bit number by combinging that last 3 of three
		uint64_t rand =
			 ((gen() & 0x0000000000000007) << 6)
			| ((gen() & 0x0000000000000007) << 3)
			| (gen() & 0x0000000000000007);
		// Only output the last 8 bits
		std::cout << reinterpret_cast<char*>(&rand)[0];
	}
	
	return 0;
}
