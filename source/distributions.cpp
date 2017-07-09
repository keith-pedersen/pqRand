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

typename pqRand::real_t pqRand::standard_normal::operator()(pqRand::engine& gen)
{
	if(valueCached)
	{
		valueCached = false;
		return cache;
	}
	else
	{
		valueCached = true;
		two const pair = GenTwo(gen);
		cache = pair.y;
		return pair.x;
	}
}

////////////////////////////////////////////////////////////////////////

typename pqRand::two pqRand::standard_normal::GenTwo(pqRand::engine& gen)
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
		pair.x = gen.U_Q();
		pair.y = gen.U_Q();
	
		// Draw x and y from U, and reject when we don't land in the circle
		u = pair.x*pair.x + pair.y*pair.y;
		
		// Reject 2/3 of the region that rounds to 1
		// (we want the epsilon/2 to the left of 1, but not the epsilon the right of 1).
		// There is a small region near (1, 0), (0, 1), etc. where the right region 
		// doesn't exist, but it is vanishingly small).
		if((u == real_t(1.)) and (gen.U_S()*real_t(3.) < real_t(2.)))
			u = 2.; // Easy way to reject
	}
	while(u > real_t(1.));

	// Give x and y a random sign via the pqRand::engine (using its bitCache)
	gen.ApplyRandomSign(pair.x);
	gen.ApplyRandomSign(pair.y);
	
	// Implement the quantile flip-flop and scale x and y
	u = (gen.RandBool() ? 
		std::sqrt(real_t(-2.)*std::log(real_t(0.5)*u)/u) : 
		std::sqrt(real_t(-2.)*std::log1p(real_t(-0.5)*u)/u));		
	
	pair.x *= u;
	pair.y *= u;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Normal takes numbers from standard_normal and adjusts them to sigma and mu
typename pqRand::two pqRand::normal::GenTwo(pqRand::engine& gen)
{
	two pair = pqRand::standard_normal::GenTwo(gen);
	
	pair.x = mu_ + sigma_ * pair.x;
	pair.y = mu_ + sigma_ * pair.y;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// Log_normal also takes from standard_normal, then exponentiates
typename pqRand::two pqRand::log_normal::GenTwo(pqRand::engine& gen)
{
	// Draw from standard normal, apply mu as a multiplicative scale
	two pair = pqRand::standard_normal::GenTwo(gen);
	
	pair.x = muScale * std::exp(sigma_ * pair.x);
	pair.y = muScale * std::exp(sigma_ * pair.y);
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// When we implement the quantile flip-flop, we must draw from HalfU
typename pqRand::real_t pqRand::weibull::operator()(pqRand::engine& gen)
{
	real_t const hu = gen.HalfU_Q();
	
	if(gen.RandBool())
		return lambda_ * std::pow(-std::log(hu), kRecip);
	else
		return lambda_ * std::pow(-std::log1p(-hu), kRecip);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::pareto::operator()(pqRand::engine& gen)
{
	return xm_ * std::pow(gen.U_Q(), negRecipAlpha);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqRand::real_t pqRand::exponential::operator()(pqRand::engine& gen)
{
	real_t const hu = gen.HalfU_Q();
	
	if(gen.RandBool())
		return -std::log(hu)/lambda_;
	else
		return -std::log1p(-hu)/lambda_;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// This version will draw 0 occasionally (when either x or y is 0, but never both).
typename pqRand::two pqRand::standard_normal_lowPrecision::GenTwo(pqRand::engine& gen)
{
	two pair;
	real_t u;
	
	do
	{
		pair.x = real_t(1.) - real_t(2.) * gen.U_S();
		pair.y = real_t(1.) - real_t(2.) * gen.U_S();
	
		// Draw x and y from U, and reject when we don't land in the circle
		u = pair.x*pair.x + pair.y*pair.y;
	}
	while(u >= real_t(1.) or (u == 0.));

	u = std::sqrt(real_t(-2.)*std::log(u)/u);
	
	pair.x *= u;
	pair.y *= u;
	
	return pair;
}
