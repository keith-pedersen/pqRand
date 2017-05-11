#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include "pqRand.hpp"
#include <cmath> // exp

namespace pqRand
{
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
			
			real_t operator()(pqRand::engine& gen); // Return a single value
			virtual two GenTwo(pqRand::engine& gen); // Return a pair of values
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
			virtual two GenTwo(pqRand::engine& gen);
			
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
				
			virtual two GenTwo(pqRand::engine& gen);
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
				lambda_(lambda_in), k_(k_in), kRecip(real_t(1.)/k_) {}
				
			real_t operator()(pqRand::engine& gen);
			
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
				xm_(xm_in), alpha_(alpha_in), negRecipAlpha(real_t(-1.)/alpha_) {}
				
			real_t operator()(pqRand::engine& gen);
			
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
				
			real_t operator()(pqRand::engine& gen);
			
			inline real_t mu() {return mu_;}
	};
	
	// Included for testing/validation purposes
	class standard_normal_lowPrecision : public standard_normal
	{
		public:
			standard_normal_lowPrecision():
				standard_normal() {}
			
			virtual two GenTwo(pqRand::engine& gen);
	};
}

#endif
