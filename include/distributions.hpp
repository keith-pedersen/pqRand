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

/*!
 *  @file distributions.hpp
 *  @brief This file Defines pqRand's non-uniform distributions.
 * 
 *  Most distributions utilize a quantile flip-flop
 *  (as described in the \ref theory "research paper")
 *  but some do not need to (e.g. pqRand::Pareto ).
 * 
 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
 *  @date 2017
*/

#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include "pqRand.hpp"
#include <cmath> // exp

namespace pqRand
{
	/*! @brief A simple struct to store a pair of doubles.
	 * 
	 *  Used by the normal distributions, since the 
	 *  Marsaglia polar generates two numbers per call.
	 * 
	 *  \warning The default constructor does not initialize its members.
	*/ 
	struct two
	{
		real_t x;
		real_t y;
		
		//! @brief Initialized construction
		two(real_t const x_in, real_t const y_in):
			x(x_in), y(y_in) {}
		
		two() {} //!< @brief \b uninitialized construction
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Samples from the standard normal distribution (\f$ \mu =0 \f$, \f$\sigma=1\f$):
	 *  \f$ \text{PDF} = \exp(-\frac{1}{2}x^2)/\sqrt{2\pi} \f$
	 * 
	 *  Samples are generated using the Marsaglia polar method, 
	 *  modified to use a quantile flip-flop and quasiuniform \f$ U \f$.
	 * 
	 *  \note The Marsaglia polar method generates two numbers per call. 
	 *  If only one number is requested, the second is cached for the next call. 
	 *  This caching functionality is the reason that this is an object, 
	 *  versus a standalone function (since its mean and variance are hard-coded). 
	 *  Furthermore, derivative distributions (\ref normal, log_normal)
	 *  extend standard_normal to reuse its caching functionality.
	*/ 
	class standard_normal
	{
		private:
			real_t cache; // Marsaglia polar samples two, cache one if we only request one
			bool valueCached; // if(valueCached == true), cache holds the next value
			
		public:
			//! @brief The distribution is hard-coded; no arguments to supply.
			standard_normal():
				valueCached(false) {}
			
			real_t operator()(pqRand::engine& gen); //!< Sample one value
			
			virtual two GenTwo(pqRand::engine& gen); //!< @brief Sample a pair of values
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Samples from the normal distribution (mean \f$ \mu \f$ and 
	 *  standard deviation \f$ \sigma \f$):
	 *  \f$ \text{PDF} = \frac{1}{\sigma\sqrt{2\pi}}
	 *  \exp\left(-\frac{1}{2\,\sigma^2}(x - \mu)^2\right) \f$
	 * 
	 *  \ref normal extends standard_normal by redefining GenTwo() to return 
	 *  \f$ x = \mu + \sigma\, x_{\mathrm{standard}}^{}  \f$.
	*/ 
	class normal : public standard_normal
	{
		protected:
			real_t const mu_;
			real_t const sigma_;
			
		public:
			normal(real_t const mu_in, real_t const sigma_in):
				standard_normal(), mu_(mu_in), sigma_(sigma_in) {}
			
			virtual two GenTwo(pqRand::engine& gen);
			
			inline real_t Mu() {return mu_;} //!< mean
			inline real_t Sigma() {return sigma_;} //!< standard deviation
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Samples from the log-normal distribution:
	 *  \f$ \text{PDF} = \frac{1}{x\,\sigma\sqrt{2\pi}}
	 *  \exp\left(-\frac{1}{2\,\sigma^2}(\log x - \mu)^2\right) \f$
	 * 
	 *  log_normal extends normal by redefining GenTwo() to return 
	 *  \f$ x = \exp (x_{\mathrm{normal}}^{}) \f$.
	*/ 
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
	
	/*! @brief Samples from the weibull distribution
	 *  (with shape parameter \f$ k \f$ and scale parameter \f$ \lambda \f$):
	 *  \f$ \text{PDF} = \frac{k}{\lambda} \left(\frac{x}{\lambda}\right)^{k-1}
	 *  \exp\left(-(x/\lambda)^k\right)\f$   
	 * 
	 *  This distribution has an invertible CDF:
	 *  \f$ \text{CDF}^{-1} = \lambda(-\log(1-u))^{1/k} \f$.
	 *  This \f$ \mathrm{CDF}^{-1} \f$ becomes ill-conditioned 
	 *  in both tails, and therefore uses a quantile flip-flop.
	*/ 
	class weibull
	{
		private:
			real_t const lambda_;
			real_t const k_;
			real_t const kRecip;
			
		public:
			weibull(real_t const lambda_in, real_t const k_in):
				lambda_(lambda_in), k_(k_in), kRecip(real_t(1.)/k_) {}
				
			real_t operator()(pqRand::engine& gen); //!< Sample one value
			
			inline real_t Lambda() {return lambda_;} //!< scale parameter
			inline real_t k() {return k_;} //!< shape parameter
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Samples from the Pareto distribution (for minimum value 
	 *  \f$ x_m^{} \f$ and index \f$ \alpha \f$): 
	 *  \f$ \text{PDF} = (\alpha\,x_m^\alpha)/(x^{\alpha  + 1}) \f$
	 * 
	 *  This distribution has an invertible CDF:
	 *  \f$ \text{CDF}^{-1} = x_m^{} u^{-1/\alpha} \f$.
	 *  This \f$ \mathrm{CDF}^{-1} \f$ becomes ill-conditioned 
	 *  in only one tail \f$(u\to 0)\f$, so we don't need a quantile flip-flop. 
	 *  However, we still need to use quasiuniform \f$ U(0, 1] \f$.
	*/ 
	class pareto
	{
		private:
			real_t const xm_;
			real_t const alpha_;
			real_t const negRecipAlpha;
			
		public:
			pareto(real_t const xm_in, real_t const alpha_in):
				xm_(xm_in), alpha_(alpha_in), negRecipAlpha(real_t(-1.)/alpha_) {}
				
			real_t operator()(pqRand::engine& gen); //!< Sample one value
			
			inline real_t xM() {return xm_;} //!< min value
			inline real_t Alpha() {return alpha_;} //!< Pareto index
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	/*! @brief Samples from the exponential distribution 
	 *  (with rate parameter \f$ \lambda \f$):
	 *  \f$ \text{PDF} = \lambda \exp(-\lambda\,x) \f$
	 * 
	 *  This distribution has an invertible CDF: 
	 *  \f$ \text{CDF}^{-1} = - \frac{1}{\lambda}\log(1-u) \f$.
	 *  This \f$ \mathrm{CDF}^{-1} \f$ becomes ill-conditioned 
	 *  in both tails, and therefore benefits from a quantile flip-flop.
	*/
	class exponential
	{
		private:
			real_t const lambda_;
			
		public:
			exponential(real_t const lambda_in):
				lambda_(lambda_in) {}
				
			real_t operator()(pqRand::engine& gen); //!< Sample one value
			
			inline real_t Lambda() {return lambda_;} //!< rate parameter
	};
	
	/*! @brief Samples from the normal distribution using the canonical polar method 
	 *  (lower precision standard_normal).
	 * 
	 *  This implementation of Marsaglia polar method is used by 
	 *  GNU and Numpy, and is included for testing purposes.
	*/
	class standard_normal_lowPrecision : public standard_normal
	{
		public:
			standard_normal_lowPrecision():
				standard_normal() {}
			
			virtual two GenTwo(pqRand::engine& gen);
	};
}

#endif
