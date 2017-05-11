#include <random>
#include <fstream>
#include <iostream>

int main()
{
	std::random_device gen;
	
	// 	sudo yum install dieharder
	// Continously output random numbers to std::out, for use like ... 
	// 	./random_device_tester.x | dieharder -g 200 -a -k 2
	while(true)
	{
		uint64_t rand = gen();
		for(size_t i = 0; i < 4; ++i)
			std::cout << reinterpret_cast<char*>(&rand)[i];
	}
	
	return 0;
}
