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

#include "../include/distributions.hpp"

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template<typename T>
std::vector<T> pqRand::distribution<T>::GetSample(size_t const sampleSize, engine& gen) const
{
	std::vector<T> sample;
	sample.reserve(sampleSize);
	
	for(size_t i = 0; i < sampleSize; ++i)
		sample.push_back((*this)(gen));

	return sample;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template<typename int_t>
pqRand::uniform_integer<int_t>::uniform_integer(int_t const min_in, int_t const max_in):
	min_(min_in), max_(max_in), 
	spread(rand_t(max_ - min_)), // The number of members of the sample space
	
	// The number of members of the PRNG sample space is (biggestRand + 1).
	// To uniformly sample, we must discard the remainder of (biggestrand + 1)/spread, 
	// so maxRand = (biggestRand - ((biggestRand + 1) % spread).
	// But (biggestRand + 1) is *potentially* not representable, 
	// so we use a trick of modular arithmetic
	// 	(a + b) % c == ((a % c) + (b % c)) % c
	maxRand(biggestRand - ((biggestRand % spread) + 1) % spread)
{
	if(max_ <= min_)
		throw std::domain_error("pqRand::uniformInteger: max must be greater than min");
	if(spread > biggestRand)
		throw std::domain_error("pqRand::uniformInteger: max is too large \
		(since we must discard the lowest engine::badBits from the PRNG)");
}

////////////////////////////////////////////////////////////////////////

template<typename int_t>
int_t pqRand::uniform_integer<int_t>::operator()(pqRand::engine& gen) const
{
	rand_t x;
	while((x = (gen() >> rightShift)) > maxRand);
	
	return int_t(x % spread) + min_;
}

// Instantiate the common types
template class pqRand::uniform_integer<int32_t>;
template class pqRand::uniform_integer<int64_t>;
template class pqRand::uniform_integer<uint32_t>;
template class pqRand::uniform_integer<uint64_t>;

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::distributionPDF::PDF(real_t const x) const
{
	if((x >= min()) and (x <= max()))
		return this->PDF_supported(x);
	else 
		return real_t(0);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::two pqRand::MeanAndVariance(distributionPDF const& dist,
	size_t const sampleSize, engine& gen)
{
	real_t sum = real_t(0);
	real_t sum2 = real_t(0);
	
	{
		auto const sample = dist.GetSample(sampleSize, gen);
						
		for(real_t const x : sample)
		{
			sum += x;
			sum2 += Squared(x);
		}
	}
	
	real_t const n = real_t(sampleSize);
	
	// This may not be the most accurate method to get the standard deviation using floats
	return two(sum / n, sum2 / n - Squared(sum) / Squared(n));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::distributionCDF::CDF(real_t const x) const
{
	if(x <= min()) return real_t(0);
	else if (x >= max()) return real_t(1);
	else return this->CDF_small_supported(x);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::distributionCDF::CompCDF(real_t const x) const
{
	if(x >= max()) return real_t(0);
	else if (x <= min()) return real_t(1);
	else return this->CDF_large_supported(x);
}
			
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::distributionQ2::operator()(engine& gen) const
{
	if(gen.RandBool())
		return Q_small(gen.HalfU_uneven());
	else
		return Q_large(gen.HalfU_uneven());
}

////////////////////////////////////////////////////////////////////////

typename pqRand::two pqRand::distributionQ2::GetTwo_antithetic(engine& gen) const
{
	real_t const hu = gen.HalfU_uneven();
	
	return two(Q_small(hu), Q_large(hu));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::uniform::uniform(real_t const min_in, real_t const max_in):
	min_(min_in), spread(max_in - min_in)
{
	if(spread <= real_t(0)) throw std::domain_error("pqRand::uniform: max must be greater than min.");
}

////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:29)
typename pqRand::real_t pqRand::uniform::PDF_supported(real_t const x) const 
{
	return real_t(1)/spread;
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::uniform::CDF_small_supported(real_t const x) const 
{
	return (x - min())/spread;
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::uniform::CDF_large_supported(real_t const x) const 
{
	return (max() - x)/spread;
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::uniform::operator()(pqRand::engine& gen) const
{
	return min_ + spread * gen.U_uneven();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:33)
pqRand::real_t pqRand::standard_normal::PDF_supported(real_t const x) const
{
	return std::exp(-real_t(0.5)*Squared(x))/real_t(std::sqrt(2 * M_PI));
}

////////////////////////////////////////////////////////////////////////
				
pqRand::real_t pqRand::standard_normal::CDF_small_supported(real_t const x) const
{
	// To use the error function we have to add, creating a cancellation
	return CDF_large_supported(-x);
}

////////////////////////////////////////////////////////////////////////

pqRand::real_t pqRand::standard_normal::CDF_large_supported(real_t const x) const
{
	return real_t(0.5)*std::erfc(x/std::sqrt(real_t(2)));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::standard_normal::operator()(pqRand::engine& gen) const
{
	if(valueCached)
	{
		valueCached = false;
		return cache;
	}
	else
	{
		valueCached = true;
		two const pair = GetTwo(gen);
		cache = pair.y;
		return pair.x;
	}
}

////////////////////////////////////////////////////////////////////////

typename pqRand::two pqRand::standard_normal::GetTwo(pqRand::engine& gen) const
{
	two pair;
	real_t u;
	
	// The Marsaglia polar method -- modeled after GNU's std::normal_distribution 
	// (bits/random.tcc, line 1925, <https://gcc.gnu.org/onlinedocs/gcc-4.8.5/libstdc++/api/a01147_source.html>)
	// on or about March 2017 -- but with the added precision of the quantile flip-flop and U_Q.
	// As such, we don't do x = (1 - U) to get a cheap random sign; 
	// this destroys the precision of quasiuniform U.
	do
	{
		pair.x = gen.U_uneven();
		pair.y = gen.U_uneven();
	
		// Draw x and y from U, and reject when we don't land in the circle
		u = pair.x*pair.x + pair.y*pair.y;
		
		// Reject 2/3 of the region that rounds to 1
		// (we want the epsilon/2 to the left of 1, but not the epsilon the right of 1).
		// There is a small region near (1, 0), (0, 1), etc. where the right region 
		// doesn't exist, but it is vanishingly small).
		if((u == real_t(1)) and (gen.U_even()*real_t(3) < real_t(2)))
			u = 2.; // Easy way to reject
	}
	while(u > real_t(1));

	// Give x and y a random sign via the pqRand::engine (using its bitCache)
	gen.ApplyRandomSign(pair.x);
	gen.ApplyRandomSign(pair.y);
	
	// Implement the quantile flip-flop and scale x and y
	u = (gen.RandBool() ? 
		std::sqrt(-real_t(2)*std::log(real_t(0.5)*u)/u) : 
		std::sqrt(-real_t(2)*std::log1p(-real_t(0.5)*u)/u));		
	
	pair.x *= u;
	pair.y *= u;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////

std::vector<pqRand::real_t> pqRand::standard_normal::GetSample
	(size_t const sampleSize, pqRand::engine& gen) const
{
	size_t const loopSize = (sampleSize + 1)/2; // This is at least sampleSize / 2;
			
	std::vector<real_t> sample;
	sample.reserve(sampleSize);
	
	for(size_t i = 0; i < loopSize; ++i)
	{
		two const pair = this->GetTwo(gen);
		sample.push_back(pair.x);
		sample.push_back(pair.y);
	}
	
	// sample.size() is now a multiple of 2; make it odd if requested
	if(sample.size() > sampleSize) 
		sample.pop_back();
	assert(sample.size() == sampleSize);
	
	return sample;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::normal::normal(real_t const mu_in, real_t const sigma_in):
	standard_normal(), mu_(mu_in), sigma_(sigma_in) 
{
	if(sigma_ <= real_t(0)) throw std::domain_error("pqRand::normal: sigma must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

// Normal takes numbers from standard_normal and adjusts them to sigma and mu
typename pqRand::two pqRand::normal::GetTwo(pqRand::engine& gen) const
{
	two pair = pqRand::standard_normal::GetTwo(gen);
	
	pair.x = mu_ + sigma_ * pair.x;
	pair.y = mu_ + sigma_ * pair.y;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:35)
typename pqRand::real_t pqRand::normal::PDF_supported(real_t const x) const
{
	return std::exp(-real_t(0.5)*Squared(x - mu_)/Squared(sigma_))/
			(sigma_*real_t(std::sqrt(2*M_PI)));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::normal::CDF_small_supported(real_t const x) const
{
	return real_t(0.5)*std::erfc((mu_-x)/(std::sqrt(real_t(2))*sigma_));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::normal::CDF_large_supported(real_t const x) const
{
	return real_t(0.5)*std::erfc((x-mu_)/(std::sqrt(real_t(2))*sigma_));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// This is the syntax to catch an exception from the super-ctor
pqRand::log_normal::log_normal(real_t const mu_in, real_t const sigma_in)
try : normal(mu_in, sigma_in), muScale(std::exp(mu_)) {}
catch(std::domain_error const& e)
{
	// Catch a domain_error from normal, to rename it (and redirect blame)
	throw std::domain_error("pqRand::log_normal: sigma must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:36)
typename pqRand::real_t pqRand::log_normal::PDF_supported(real_t const x) const
{
	return normal::PDF_supported(std::log(x))/x;
}

////////////////////////////////////////////////////////////////////////
		
typename pqRand::real_t pqRand::log_normal::CDF_small_supported(real_t const x) const
{
	return normal::CDF_small_supported(std::log(x));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::log_normal::CDF_large_supported(real_t const x) const
{
	return normal::CDF_large_supported(std::log(x));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::log_normal::Mean() const 
{
	return std::exp(mu_ + real_t(0.5)*Squared(sigma_));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::log_normal::Variance() const 
{
	return std::exp(real_t(2)*mu_ + Squared(sigma_))*std::expm1(Squared(sigma_));
}

////////////////////////////////////////////////////////////////////////

// Log_normal also takes from standard_normal, then exponentiates
typename pqRand::two pqRand::log_normal::GetTwo(pqRand::engine& gen) const
{
	// Draw from standard normal, apply mu as a multiplicative scale
	two pair = pqRand::standard_normal::GetTwo(gen);
	
	pair.x = muScale * std::exp(sigma_ * pair.x);
	pair.y = muScale * std::exp(sigma_ * pair.y);
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::weibull::weibull(real_t const lambda_in, real_t const k_in):
	lambda_(lambda_in), k_(k_in), kRecip(real_t(1)/k_)
{
	if(lambda_ <= real_t(0)) 
		throw std::domain_error("pqRand::weibull: lambda must be greater than zero!");
	if(kRecip <= real_t(0))
		throw std::domain_error("pqRand::weibull: k must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:39)
typename pqRand::real_t pqRand::weibull::PDF_supported(real_t const x) const
{
	real_t const xOverLambda = x / lambda_;
	real_t const xl2kminus1 = std::pow(xOverLambda, k_ - real_t(1));
	
	return (std::exp(-xl2kminus1 * xOverLambda) * xl2kminus1 * k_)/lambda_;
}

////////////////////////////////////////////////////////////////////////
					
typename pqRand::real_t pqRand::weibull::CDF_small_supported(real_t const x) const
{
	return -std::expm1(-std::pow(x/lambda_, k_));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::weibull::CDF_large_supported(real_t const x) const 
{
	return std::exp(-std::pow(x/lambda_, k_));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::weibull::Mean() const 
{
	return lambda_ * std::tgamma(real_t(1) + real_t(1) / k_);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::weibull::Variance() const 
{
	return Squared(lambda_) * 
				(std::tgamma(real_t(1) + real_t(2) / k_) - 
					Squared(std::tgamma(real_t(1) + real_t(1) / k_)));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::weibull::Q_small(real_t const u) const
{
	return lambda_ * std::pow(-std::log1p(-u), kRecip);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::weibull::Q_large(real_t const u) const
{
	return lambda_ * std::pow(-std::log(u), kRecip);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::pareto::pareto(real_t const xMin_in, real_t const alpha_in):
	xMin(xMin_in), alpha_(alpha_in), negRecipAlpha(-real_t(1)/alpha_), 
	alpha_xM2alpha(alpha_*std::pow(xMin, alpha_))
{
	if(xMin <= real_t(0)) 
		throw std::domain_error("pqRand::pareto: x_m must be greater than zero!");
	if(negRecipAlpha >= real_t(0)) 
		throw std::domain_error("pqRand::pareto: alpha must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::pareto::Mean() const 
{
	return (alpha_ <= real_t(1)) ? INFINITY : alpha_ * xMin / (alpha_ - real_t(1));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::pareto::Variance() const 
{
	return (alpha_ <= real_t(2)) ? INFINITY : 
		alpha_ * Squared(xMin) / (Squared(alpha_ - real_t(1)) * (alpha_ - real_t(2)));
}			

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::pareto::operator()(pqRand::engine& gen) const
{
	return xMin * std::pow(gen.U_uneven(), negRecipAlpha);
}

////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:43)
typename pqRand::real_t pqRand::pareto::PDF_supported(real_t const x) const
{
	return std::pow(x, -(alpha_ + real_t(1)))*alpha_xM2alpha;
}

////////////////////////////////////////////////////////////////////////
					
typename pqRand::real_t pqRand::pareto::CDF_small_supported(real_t const x) const 
{
	return -std::expm1(alpha_ * std::log(xMin/x));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::pareto::CDF_large_supported(real_t const x) const
{
	return std::exp(alpha_ * std::log(xMin/x));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::exponential::exponential(real_t const lambda_in):
	lambda_(lambda_in)
{
	if(lambda_ <= real_t(0))
		throw std::domain_error("pqRand::exponential: lambda must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:43)
typename pqRand::real_t pqRand::exponential::PDF_supported(real_t const x) const
{
	return lambda_ * std::exp(-lambda_ * x);
}

////////////////////////////////////////////////////////////////////////
					
typename pqRand::real_t pqRand::exponential::CDF_small_supported(real_t const x) const
{
	return -std::expm1(-lambda_ * x);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::exponential::CDF_large_supported(real_t const x) const
{
	return std::exp(-lambda_ * x);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::exponential::Q_small(real_t const u) const
{
	return -std::log1p(-u)/lambda_;
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::exponential::Q_large(real_t const u) const
{
	return -std::log(u)/lambda_;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::logistic::logistic(real_t const mu_in, real_t const s_in):
	mu_(mu_in), s_(s_in)
{
	if(s_ <= real_t(0))
		throw std::domain_error("pqRand::exponential: s must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

// Checked PDF/CDF (16.12.2017 @ 11:44)
typename pqRand::real_t pqRand::logistic::PDF_supported(real_t const x) const
{
	real_t const expTerm = std::exp(-(x-mu_)/s_);
	
	return expTerm/(s_ * Squared(real_t(1) + expTerm));
}

////////////////////////////////////////////////////////////////////////
								
typename pqRand::real_t pqRand::logistic::CDF_small_supported(real_t const x) const
{
	return real_t(1) / (real_t(1) + std::exp(-(x-mu_)/s_));
}

////////////////////////////////////////////////////////////////////////
	
typename pqRand::real_t pqRand::logistic::CDF_large_supported(real_t const x) const
{
	return real_t(1) / (real_t(1) + std::exp((x-mu_)/s_));
}

////////////////////////////////////////////////////////////////////////

//~ return mu_ + s_ * gen.ApplyRandomSign(std::log(real_t(1)/gen.HalfU_uneven() - real_t(1)));

typename pqRand::real_t pqRand::logistic::Q_small(real_t const u) const
{
	return mu_ - s_ * std::log(real_t(1)/u - real_t(1));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::logistic::Q_large(real_t const u) const
{
	return mu_ + s_ * std::log(real_t(1)/u - real_t(1));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::log_logistic::log_logistic(real_t const alpha_in, real_t const beta_in):
	alpha_(alpha_in), beta_(beta_in), betaInverse(real_t(1)/beta_)
{
	if(alpha_ <= real_t(0))
		throw std::domain_error("pqRand::exponential: alpha must be greater than zero!");
	if(beta_ <= real_t(0))
		throw std::domain_error("pqRand::exponential: beta must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::log_logistic::PDF_supported(real_t const x) const
{
	real_t const xOverAlpha = x/alpha_;
	real_t const xOverAlpha2betaMinus1 = std::pow(xOverAlpha, beta_ - real_t(1));
					
	return (xOverAlpha2betaMinus1 * beta_) / 
		(alpha_ * Squared(real_t(1) + xOverAlpha2betaMinus1 * xOverAlpha));
}

////////////////////////////////////////////////////////////////////////
					
typename pqRand::real_t pqRand::log_logistic::CDF_small_supported(real_t const x) const
{
	return real_t(1) / (real_t(1) + std::pow(x/alpha_, -beta_));
}

////////////////////////////////////////////////////////////////////////
	
typename pqRand::real_t pqRand::log_logistic::CDF_large_supported(real_t const x) const
{
	return real_t(1) / (real_t(1) + std::pow(x/alpha_, beta_));
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::log_logistic::Mean() const 
{
	return (beta_ <= real_t(1)) ? INFINITY : 
	(alpha_*M_PI)/(beta_ * std::sin(M_PI / beta_));
}

////////////////////////////////////////////////////////////////////////
		
typename pqRand::real_t pqRand::log_logistic::Variance() const
{
	if (beta_ <= real_t(2)) return INFINITY;
	
	real_t const b = M_PI / beta_;
	
	return Squared(alpha_) * (b * (real_t(1)/std::cos(b) - b / std::sin(b)))/std::sin(b);
}

////////////////////////////////////////////////////////////////////////

//~ return alpha_ * std::pow(real_t(1)/gen.HalfU_uneven() - real_t(1), 
		//~ gen.ApplyRandomSign(real_t(betaInverse))); // real_t to make a copy we can alter
		
typename pqRand::real_t pqRand::log_logistic::Q_small(real_t const u) const
{
	return alpha_ * std::pow(real_t(1)/u - real_t(1), -betaInverse);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::log_logistic::Q_large(real_t const u) const
{
	return alpha_ * std::pow(real_t(1)/u - real_t(1), betaInverse);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

pqRand::gammaDist::gammaDist(real_t const lambda_in, real_t const k_in)
	try : 
	lambda_(lambda_in), k_(k_in), 
	proposal(k_, std::sqrt(real_t(2)*k_ - real_t(1))),
	lambda2k(std::pow(lambda_, k_)), logGamma_k(std::lgamma(k_))
{
	if(lambda_ <= real_t(0))
		throw std::domain_error("pqRand::gammaDist: lambda must be greater than zero!");
	if(k_ <= real_t(1))
		throw std::domain_error("pqRand::gammaDist: k must be greater than 1!");
}
catch(std::domain_error const& e)
{
	// log_logistic will catch this first
	throw std::domain_error("pqRand::gammaDist: k must be greater than zero!");
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::gammaDist::PDF_supported(real_t const x) const
{
	return lambda2k * 
		std::exp((k_ - real_t(1)) * std::log(x) - lambda_ * x - logGamma_k);
}

////////////////////////////////////////////////////////////////////////

bool pqRand::gammaDist::Reject(real_t const x, pqRand::engine& gen) const
{
	real_t const xOverK = x / k_;
	real_t const xTerm = std::pow(xOverK, proposal.Beta());
	
	// The simplified ratio of the target PDF to the proposal PDF
	real_t const acceptProbability = 
	((real_t(0.25)*std::exp(k_ - x)*std::pow(xOverK, k_)*
		Squared(real_t(1) + xTerm)) / xTerm);
	
	assert((acceptProbability - real_t(1)) < real_t(1e-8));
	
	// If a uniform variate is greater than the acceptance rate, then reject
	return (gen.U_uneven() > acceptProbability);
}

////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::gammaDist::operator()(pqRand::engine& gen) const
{
	// sample from lambda = 1, then scale by the actual lambda.
	real_t x;
	while(Reject(x = proposal(gen), gen));
	
	return x / lambda_;
}
			
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// This version will draw 0 occasionally (when either x or y is 0, but never both).
typename pqRand::two pqRand::standard_normal_lowPrecision::GetTwo(pqRand::engine& gen) const
{
	two pair;
	real_t u;
	
	do
	{
		pair.x = real_t(1.) - real_t(2.) * gen.U_even();
		pair.y = real_t(1.) - real_t(2.) * gen.U_even();
	
		// Draw x and y from U, and reject when we don't land in the circle
		u = pair.x*pair.x + pair.y*pair.y;
	}
	while(u >= real_t(1.) or (u == 0.));

	u = std::sqrt(real_t(-2.)*std::log(u)/u);
	
	pair.x *= u;
	pair.y *= u;
	
	return pair;
}
