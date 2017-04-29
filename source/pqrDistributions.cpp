#include "pqrDistributions.hpp"

typename pqr::real_t pqr::standard_normal::operator()(pqRand_engine& gen)
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
		cache = pair.x;
		return pair.y;
	}
}

////////////////////////////////////////////////////////////////////////

typename pqr::two pqr::standard_normal::GenTwo(pqRand_engine& gen)
{
	two pair;
	real_t u;
	
	// The Marsaglia polar method. Note that we're not doing 1 - U
	// to get a cheap random sign; this destroys the precision of U
	do
	{
		pair.x = gen.U_Q();
		pair.y = gen.U_Q();
	
		// Draw x and y from U, and reject when we don't land in the circle
		u = pair.x*pair.x + pair.y*pair.y;
		
		// Reject 2/3 of the region that rounds to 1
		// (we want the epsilon/2 to the left of 1, but not the epsilon the right of 1.
		//  There is a small region near (1, 0), (0, 1), etc. where the R region 
		//  doesn't exist, but it is vanishingly small).
		if((u == 1.) and (gen.U_Q()*3. < 2.))
			u = 2.; // Easy way to reject
	}
	while(u > 1.);

	// Give x and y a random sign via the pqRand_engine (uses its bitCache)
	gen.ApplyRandomSign(pair.x);
	gen.ApplyRandomSign(pair.y);
	
	// Implement the quantile flip-flop and scale x and y
	u = (gen.RandBool() ? 
		std::sqrt(-2.*std::log(0.5*u)/u) : 
		std::sqrt(2.*std::log1p(u/(2. - u))/u));		
	
	pair.x *= u;
	pair.y *= u;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Normal takes numbers from standard_normal and adjusts them to sigma and mu
typename pqr::two pqr::normal::GenTwo(pqRand_engine& gen)
{
	two pair = pqr::standard_normal::GenTwo(gen);
	
	pair.x = mu_ + sigma_ * pair.x;
	pair.y = mu_ + sigma_ * pair.y;
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// Log_normal also takes from standard_normal, then exponentiates
typename pqr::two pqr::log_normal::GenTwo(pqRand_engine& gen)
{
	// Draw from standard normal, apply mu as a multiplicative scale
	two pair = pqr::standard_normal::GenTwo(gen);
	
	pair.x = muScale * std::exp(sigma_ * pair.x);
	pair.y = muScale * std::exp(sigma_ * pair.y);
	
	return pair;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// When we implement the quantile flip-flop, we must draw from HalfU
typename pqr::real_t pqr::weibull::operator()(pqRand_engine& gen)
{
	real_t const hu = gen.HalfU_Q();
	
	if(gen.RandBool())
		return lambda_ * std::pow(-std::log(hu), kRecip);
	else
		return lambda_ * std::pow(std::log1p(hu/(1.-hu)), kRecip);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqr::real_t pqr::pareto::operator()(pqRand_engine& gen)
{
	return xm_ * std::pow(gen.U_Q(), negRecipAlpha);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typename pqr::real_t pqr::exponential::operator()(pqRand_engine& gen)
{
	real_t const hu = gen.HalfU_Q();
	
	if(gen.RandBool())
		return -mu_ * std::log(hu);
	else
		return mu_ * std::log1p(hu/(1.-hu));
}
