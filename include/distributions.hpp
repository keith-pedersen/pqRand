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
*  @brief Defines the sampling objects which draw from various distributions.
* 
*  The sampling objects use a polymorphic interface which
*  adds functions as more analytic information is known about their distribution.
*  Most distributions utilize a quantile flip-flop (pqRand::distributionQ2,
*  as described in the \ref theory "research paper")
*  but some do not need one (e.g. pqRand::pareto).
*  Others use rejection sampling (e.g. pqRand::standard_normal, pqRand::gammaDist).
* 
*  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
*  @date 2017
*/

#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include "pqRand.hpp"
#include <cmath> // exp
#include <assert.h>

namespace pqRand
{
	/*! @brief An abtract (pure virtual) class for sampling random variates.
	 * 
	 *  This class defines a standard interface: the operator() is used to sample one variate, 
	 *  of which the min/max variate must be defined. 
	 *  This permits a DRY (don't repeat yourself) GetSample().
	*/
	template<typename T>
	class distribution
	{
		public:
			distribution() {}
			virtual ~distribution() {}
			
			virtual T min() const = 0; //!< @brief The minimum variate sampled.
			virtual T max() const = 0; //!< @brief The maximum variate sampled.
			
			/*! @brief Sample one variate.
			 *  
			 *  \param gen 	the PRNG engine
			*/ 
			virtual T operator()(engine& gen) const = 0;
			
			/*! @brief Sample a number of variates and return them in a vector.
			 *  
			 *  \param sampleSize 	the sample size
			 *  \param gen 	the PRNG engine
			*/ 
			virtual std::vector<T> GetSample(size_t const sampleSize, engine& gen) const;
	};	
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief A simple struct to store a pair of real_t.
	 * 
	 *  Used by std_normal and its children, since the 
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
		
		two() {} //!< @brief \b Uninitialized construction
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief An abstract class which extends the \ref distribution interface for 
	 *  those distributions with an analytic PDF.
	 * 
	 *  This allows access to the \ref PDF, the \ref Mean, and \ref Variance.
	*/
	class distributionPDF : public distribution<real_t>
	{
		protected:
			//! @brief The PDF, where x is guaranteed to be supported by \ref PDF.
			virtual real_t PDF_supported(real_t const x) const = 0;
					
		public:
			distributionPDF() {}
			virtual ~distributionPDF() {}
			
			//! @brief The probability distribution function (zero outside of [min, max]).
			real_t PDF(real_t const x) const;
			
			//! @brief The distribution's mean, \f$ \langle x \rangle = \int \text{PDF}(x)\, x\, dx \f$.
			virtual real_t Mean() const = 0; 
			//! @brief The distribution's variance, \f$ \langle x^2 \rangle - \langle x \rangle^2 \f$.
			virtual real_t Variance() const = 0;
			
			using distribution<real_t>::GetSample; // Declare "using" to force creation of binary code
	};
	
	/*! @brief Sample many variates and calculate their mean and variance (for validation).
	 * 
	 * \param sampleSize 	the sample size
	 * \param gen 		the PRNG engine
	 *
	 * \return the mean and variance in a \ref two (two.x = mean, two.y = variance).
	*/			
	two MeanAndVariance(distributionPDF const& dist, size_t const sampleSize, engine& gen);

	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief An abstract class which extends the \ref distributionPDF interface 
	 *  for those distributions with a calculable CDF.
	*/ 
	class distributionCDF : public distributionPDF
	{
		protected:
			//! @brief The CDF, where x is guaranteed to be supported by \ref CDF.
			virtual real_t CDF_small_supported(real_t const x) const = 0;
			//! @brief The complementary CDF, where x is guaranteed to be supported by \ref CDF.
			virtual real_t CDF_large_supported(real_t const x) const = 0;
		
		public:
			distributionCDF() {}		
			virtual ~distributionCDF() {}
			
			//! @brief The cumulative distribution function; \f$ \text{CDF}(x) = \int_\text{min}^x \text{PDF}(x^\prime)\text{d}x^\prime \f$.
			real_t CDF(real_t const x) const;
			
			//! @brief The accurate complementary CDF (i.e. 1 - CDF(x), without cancellation).
			real_t CompCDF(real_t const x) const;
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////	
	
	/*! @brief An abstract class which extends the \ref distributionCDF interface 
	 *  for those distributions with an invertible CDF (i.e. quantile functions).
	 * 
	 *  This permits a generic operator() using a quantile flip-flop 
	 *  (the precise inversion method, as discuseed in the \ref theory paper),
	 *  which keeps the quantile function well-conditioned in both tails.
	*/
	class distributionQ2 : public distributionCDF
	{
		public:
			distributionQ2() {}		
			virtual ~distributionQ2() {}
			
			/*! @brief The quantile function (the inverse of the CDF),
			 *  which accurately samples the small-value tail (given \p u < 1/2).
			 * 
			 *  \param u 	a uniform variate from (0, 1)
			*/ 
			virtual real_t Q_small(real_t const u) const = 0;
			
			/*! @brief The complementary quantile function 
			 *  (the inverse of the Complementary CDF, by taking u -> 1 - u in Q_small),
			 *  which accurately samples the large-value tail (given \p u < 1/2).
			 * 
			 *  \param u 	a uniform variate from (0, 1)
			*/
			virtual real_t Q_large(real_t const u) const = 0;
			
			/*! @brief Sample one variate using a quantile flip-flop.
			 * 
			 *  Randomly choose Q_small or Q_large, then feed it a variate from pqRand:engine::HalfU_uneven.
			 * 
			 * \param gen 		the PRNG engine
			*/
			real_t operator()(engine& gen) const;
			
			/*! @brief Sample antithetic variates (same u through both Q).
			 * 
			 * \param gen 		the PRNG engine
			*/	
			two GetTwo_antithetic(engine& gen) const;
	};	
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Sample integers uniformly from the half-open [min, max).
	 * 
	 *  \note max() is not actually returned; the largest sampled variate is (max - 1).
	 *  This is a slight perversion of the \ref distribution interface, 
	 *  but is most consistent with the standard definition of uniform integer sampling.
	*/	
	template<typename int_t>
	class uniform_integer : public distribution<int_t>
	{
		private:
			using rand_t = pqRand::engine::result_type;
		
			int_t const min_;
			int_t const max_;
			rand_t const spread;
			rand_t const maxRand;
			
			rand_t static constexpr rightShift = pqRand::engine::badBits;
			rand_t static constexpr biggestRand = (pqRand::engine::max() >> rightShift);
			
			// This class is simple, so it cannot create more entropy than the PRNG.
			static_assert(std::numeric_limits<int_t>::digits <= std::numeric_limits<rand_t>::digits,
				"pqRand::uniformInteger: integer type has too many digits");
			static_assert(pqRand::engine::min() == 0, 
				"pqRand::uniformInteger: PRNG must return 0.");
					
		public:
			/*! @brief Define the half-open interval [min, max)
			 *  
			 *  \throws throws std::domain_error if (\p max <= \p min).
			*/ 
			explicit uniform_integer(int_t const min, int_t const max);
			
			int_t operator()(pqRand::engine& gen) const;
			using distribution<int_t>::GetSample; // Declare "using" to force creation of binary code
			
			inline int_t min() const {return min_;}
			inline int_t max() const {return max_;} //!< One past the maximum variate sampled.
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
			
	/*! @brief Sample the uniform distribution over the \em closed interval [min, max].
	 * 
	 * \f$ \text{PDF}(x) = \frac{1}{\text{max}-\text{min}} \f$
	 * 
	 * \warning \em Both endpoints are returned (unless \p a is zero, 
	 * in which case <em> exactly zero </em> will probably not be returned).
	*/
	class uniform : public distributionCDF
	{
		private:
			real_t const min_;
			real_t const spread;
			
		protected:
			real_t PDF_supported(real_t const x) const;
			real_t CDF_small_supported(real_t const x) const;
			real_t CDF_large_supported(real_t const x) const;			
					
		public:
			/*! @brief Define the closed interval [min, max]
			 *  
			 *  \throws throws std::domain_error if (\p max <= \p min).
			*/ 
			uniform(real_t const min, real_t const max);
			
			real_t Mean() const {return real_t(0.5)*(min() + max());}
			real_t Variance() const {return Squared(max() - min())/real_t(12);}
							
			real_t operator()(pqRand::engine& gen) const;
			
			inline real_t min() const {return min_;}
			inline real_t max() const {return min_ + spread;}			
			// Note: ((max - min) + min == max) with floats if max > min. This is important.
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Sample the standard normal distribution (\f$ \mu =0 \f$, \f$\sigma=1\f$).
	 * 
	 *  \f$ \text{PDF}(x) = \exp(-\frac{1}{2}x^2)/\sqrt{2\pi} \f$
	 * 
	 *  Samples are generated using the Marsaglia polar method, 
	 *  modified to use a quantile flip-flop and uneven \f$ U(0,1) \f$.
	 * 
	 *  \note The Marsaglia polar method generates two numbers per call. 
	 *  If only one number is requested, the second is cached for the next call.
	*/ 
	class standard_normal : public distributionCDF
	{
		private:
			// Make these mutable because the object isn't really changing
			mutable real_t cache; // Marsaglia polar samples two, cache one if we only request one
			mutable bool valueCached; // if(valueCached == true), cache holds the next variate
			
		protected:
			virtual real_t PDF_supported(real_t const x) const;							
			virtual real_t CDF_small_supported(real_t const x) const;
			virtual real_t CDF_large_supported(real_t const x) const;
			
		public:
			// The distribution is hard-coded; no arguments to supply.
			standard_normal(): valueCached(false) {}
			virtual ~standard_normal() {}
			
			virtual inline real_t min() const {return -INFINITY;}
			virtual inline real_t max() const {return INFINITY;}
			
			real_t Mean() const {return real_t(0);}
			real_t Variance() const {return real_t(1);}
			
			real_t operator()(pqRand::engine& gen) const;
			virtual two GetTwo(pqRand::engine& gen) const; //!< @brief Sample a pair of variates.
			
			// Redefine GetSample to draw two variates at a time, skipping the caching mechanism.
			std::vector<real_t> GetSample(size_t const sampleSize, pqRand::engine& gen) const;
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Sample the normal distribution (mean \f$ \mu \f$ and 
	 *  standard deviation \f$ \sigma > 0 \f$).
	 * 
	 *  \f$ \text{PDF}(x) = \frac{1}{\sigma\sqrt{2\pi}}
	 *  \exp\left(-\frac{1}{2\,\sigma^2}(x - \mu)^2\right) \f$
	 * 
	 *  \ref normal extends standard_normal by redefining GetTwo() to return 
	 *  \f$ x = \mu + \sigma\, x_{\mathrm{standard}}^{}  \f$.
	*/ 
	class normal : public standard_normal
	{
		protected:
			real_t const mu_;
			real_t const sigma_;
			
			virtual real_t PDF_supported(real_t const x) const;
			virtual real_t CDF_small_supported(real_t const x) const;
			virtual real_t CDF_large_supported(real_t const x) const;
						
		public:
			/*! @brief Define the distribution's parameters.
			 * 
			 * \param mu 	the mean
			 * \param sigma the standard deviation
			 * 
			 * \throws throws std::domain_error if (\p sigma <= 0).
			*/ 
			normal(real_t const mu, real_t const sigma);			
			virtual ~normal() {}
			
			real_t Mean() const {return mu_;}
			real_t Variance() const {return Squared(sigma_);}
			
			virtual two GetTwo(pqRand::engine& gen) const;
			
			virtual inline real_t Mu() const {return mu_;} //!< The mean
			virtual inline real_t Sigma() const {return sigma_;} //!< The standard deviation
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Sample the log-normal distribution (the \em log of the variates
	 *  are normally distributed, with mean \f$ \mu \f$ and standard deviation \f$ \sigma > 0 \f$).
	 * 
	 *  \f$ \text{PDF}(x) = \frac{1}{x\,\sigma\sqrt{2\pi}}
	 *  \exp\left(-\frac{1}{2\,\sigma^2}(\log x - \mu)^2\right) \f$
	 * 
	 *  log_normal extends normal by redefining GetTwo() to return 
	 *  \f$ x = \exp (x_{\mathrm{normal}}^{}) \f$.
	*/ 
	class log_normal : public normal
	{
		private:
			real_t const muScale; // Exponentiate mu once, multiply times all returned values
			// This is done to prevent assumed addition cancellation in 
			// exp(mu + sigma * x), but is the cancellation actually avoided?
			
		protected:
			real_t PDF_supported(real_t const x) const;
			real_t CDF_small_supported(real_t const x) const;
			real_t CDF_large_supported(real_t const x) const;
			
		public:
			/*! @brief Define the distribution's parameters.
			 * 
			 * \param mu 	the mean of the variates' log
			 * \param sigma 	the standard deviation of the variates' log
			 * 
			 * \throws throws std::domain_error if (\p sigma <= 0).
			*/ 
			log_normal(real_t const mu, real_t const sigma);
			
			inline real_t min() const {return 0;}
			
			real_t Mean() const; 
			real_t Variance() const;
			
			two GetTwo(pqRand::engine& gen) const;
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Sample the Weibull distribution
	 *  (with scale \f$ \lambda > 0 \f$ and shape \f$ k > 0 \f$).
	 * 
	 *  \f$ \text{PDF}(x) = \frac{k}{\lambda} \left(\frac{x}{\lambda}\right)^{k-1}
	 *  \exp\left(-(x/\lambda)^k\right)\f$
	*/ 
	class weibull : public distributionQ2
	{
		private:
			real_t const lambda_;
			real_t const k_;
			real_t const kRecip;			
		
		protected:
			// Checked PDF/CDF (16.12.2017 @ 11:39)
			real_t PDF_supported(real_t const x) const;								
			real_t CDF_small_supported(real_t const x) const;
			real_t CDF_large_supported(real_t const x) const;
			
		public:
			/*! @brief Define the distribution's parameters.
			 * 
			 *  \param lambda 	the scale
			 *  \param k 	the shape
			 * 
			 *  \throws throws std::domain_error if either parameter is non-positive.
			*/ 
			weibull(real_t const lambda, real_t const k);
			
			inline real_t min() const {return 0;}
			inline real_t max() const {return INFINITY;}
			
			real_t Mean() const;
			real_t Variance() const;
			
			inline real_t Lambda() const {return lambda_;} //!< The scale parameter.
			inline real_t k() const {return k_;} //!< The shape parameter.
			
			real_t Q_small(real_t const u) const;
			real_t Q_large(real_t const u) const;
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	
	/*! @brief Sample the Pareto distribution (with minimum value 
	 *  \f$ x_\text{min}^{} > 0 \f$ and index \f$ \alpha > 0 \f$).
	 * 
	 *  \f$ \text{PDF}(x) = (\alpha\,x_\text{min}^\alpha)/(x^{\alpha  + 1}) \f$
	 * 
	 *  \note This distribution has an invertible CDF, but it only becomes ill-conditioned 
	 *  in one tail \f$(u\to 0)\f$, so we don't need a quantile flip-flop.
	 *  However, we still need to use uneven \f$ U(0, 1] \f$.
	*/ 
	class pareto : public distributionCDF
	{
		private:
			real_t const xMin;
			real_t const alpha_;
			real_t const negRecipAlpha;
			real_t const alpha_xM2alpha;
						
		protected:
			real_t PDF_supported(real_t const x) const;
			real_t CDF_small_supported(real_t const x) const;
			real_t CDF_large_supported(real_t const x) const;
			
		public:
			/*! @brief Define the distribution's parameters.
			 * 
			 *  \param xMin 	the minimum value
			 *  \param k 	the Pareto index
			 * 
			 *  \throws throws std::domain_error if either parameter is non-positive.
			*/ 
			pareto(real_t const xMin, real_t const alpha);
			
			inline real_t min() const {return xMin;}
			inline real_t max() const {return INFINITY;}
			
			real_t Mean() const;
			real_t Variance() const;
			
			real_t operator()(pqRand::engine& gen) const;
			
			inline real_t Alpha() const {return alpha_;} //!< The Pareto index
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	/*! @brief Sample the exponential distribution 
	 *  (with rate parameter \f$ \lambda > 0 \f$).
	 * 
	 *  \f$ \text{PDF}(x) = \lambda \exp(-\lambda\,x) \f$
	*/
	class exponential : public distributionQ2
	{
		private:
			real_t const lambda_;
			
		protected:
			real_t PDF_supported(real_t const x) const;
			real_t CDF_small_supported(real_t const x) const;
			real_t CDF_large_supported(real_t const x) const;
			
		public:
			/*! @brief Define the distribution's rate parameter \p lambda.
			 *  
			 *  \throws throws std::domain_error if (\p lambda <= 0).
			*/ 
			exponential(real_t const lambda);
			
			inline real_t min() const {return 0;}
			inline real_t max() const {return INFINITY;}
			
			real_t Mean() const {return real_t(1)/lambda_;}
			real_t Variance() const {return real_t(1)/Squared(lambda_);}
							
			real_t Q_small(real_t const u) const;
			real_t Q_large(real_t const u) const;
			
			inline real_t Lambda() const {return lambda_;} //!< The rate parameter			
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	/*! @brief Sample the logistic distribution 
	 *  (with mean \f$ \mu \f$ and scale \f$ s > 0 \f$).
	 * 
	 *  \f$ \text{PDF}(x) = \left(4s\,\cosh^2\left(\frac{x-\mu}{2s}\right)\right)^{-1} \f$
	*/
	class logistic : public distributionQ2
	{
		private:
			real_t const mu_;
			real_t const s_;
			
		protected:
			// Checked PDF/CDF (16.12.2017 @ 11:44)
			real_t PDF_supported(real_t const x) const;
			real_t CDF_small_supported(real_t const x) const;
			real_t CDF_large_supported(real_t const x) const;
			
		public:
			/*! @brief Define the distribution's parameters.
			 * 
			 *  \param mu 	the mean
			 *  \param s 	the scale
			 *  
			 *  \throws throws std::domain_error if the (\p scale <= 0).
			*/ 
			logistic(real_t const mu, real_t const s);
			
			inline real_t min() const {return -INFINITY;}
			inline real_t max() const {return INFINITY;}
			
			real_t Mean() const {return mu_;}
			real_t Variance() const {return Squared(s_ * M_PI)/real_t(3);}
							
			real_t Q_small(real_t const u) const;
			real_t Q_large(real_t const u) const;
			
			inline real_t Mu() const {return mu_;} //!< The mean
			inline real_t s() const {return s_;} //!< The scale
	};
	
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	/*! @brief Sample the log-logistic distribution 
	 *  (with scale \f$ \alpha > 0 \f$ and shape \f$ \beta > 0\f$).
	 * 
	 *  \f$ \text{PDF}(x) = \frac{(\beta/\alpha)(x/\alpha)^{\beta-1}}{(1+(x/\alpha)^\beta)^2} \f$
	*/
	class log_logistic : public distributionQ2
	{
		private:
			real_t const alpha_;
			real_t const beta_;
			real_t const betaInverse;
			
		protected:
			real_t PDF_supported(real_t const x) const;
			real_t CDF_small_supported(real_t const x) const;				
			real_t CDF_large_supported(real_t const x) const;
			
		public:
			/*! @brief Define the distribution's parameters.
			 * 
			 *  \param alpha 	the scale
			 *  \param beta 	the shape
			 *  
			 *  \throws throws std::domain_error if either parameter is non-positive.
			*/ 
			log_logistic(real_t const alpha, real_t const beta);
			
			inline real_t min() const {return 0;}
			inline real_t max() const {return INFINITY;}
			
			real_t Mean() const;
			real_t Variance() const;
			
			real_t Q_small(real_t const u) const;
			real_t Q_large(real_t const u) const;
							
			inline real_t Alpha() const {return alpha_;} //!< The scale
			inline real_t Beta() const {return beta_;} //!< The shape
	};
	
	/*! @brief Sample the gamma distribution 
	 *  (with rate \f$ \lambda > 0 \f$ and shape \f$ k > 1 \f$), 
	 *  the sum of \f$ k \f$ \ref exponential distributions with rate \f$ \lambda \f$.
	 * 
	 *  \f$ \text{PDF}(x) = \frac{\lambda^k}{\Gamma(k)}x^{k-1}\exp(-\lambda\,x) \f$
	 * 
	 *  This class uses the rejection sampling scheme proposed in
	 *  > Cheng, R.C.H. "The generation of gamma variables with non-integral shape parameter,"
	 *  > Journal of the Royal Statistical Socieity C, Vol. 26 (1977)
	 *  which samples a \f$ \lambda = 1 \f$ gamma distribution
	 *  using a special log_logistic as a proposal distribution.
	 *  The \f$ \lambda = 1 \f$ variates are internally scaled to 
	 *  any \f$ \lambda \f$ by redefining the timescale (if \f$ \lambda = 10 \f$, 
	 *  the \f$ \lambda = 1 \f$ variates are downscaled by 10).
	 * 
	 *  \note This class is named gammaDist (instead of simply \p gamma)
	 *  to prevent confusion with the gamma function.
	*/
	class gammaDist : public distributionPDF
	{
		private:
			real_t const lambda_;
			real_t const k_;
			log_logistic const proposal;
			real_t const lambda2k;
			real_t const logGamma_k;
			
			// Should the x from the proposal distribution be rejected?
			bool Reject(real_t const x, pqRand::engine& gen) const;
			
		protected:
			real_t PDF_supported(real_t const x) const;
			
		public:
			/*! @brief Define the distribution's parameters.
			 * 
			 *  \param lambda 	the rate
			 *  \param k 	the shape
			 *  
			 *  \throws throws std::domain_error if (\p lambda <= 0) or (\p k <= 1)
			*/ 
			gammaDist(real_t const lambda, real_t const k);
			
			inline real_t min() const {return 0;}
			inline real_t max() const {return INFINITY;}
			
			real_t Mean() const {return k_ / lambda_;}
			real_t Variance() const {return k_ / Squared(lambda_);}
			
			real_t operator()(pqRand::engine& gen) const;
			
			inline real_t Lambda() const {return lambda_;} //!< The rate
			inline real_t k() const {return k_;} //!< The shape
	};	
	
	/*! @brief Sample from the normal distribution using the canonical polar method 
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
				
			virtual ~standard_normal_lowPrecision() {}
			
			virtual two GetTwo(pqRand::engine& gen) const;
	};
}

#endif
