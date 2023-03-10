//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_randn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/rand/mc_randg.h>

#ifndef MC_RANDN_H
#define MC_RANDN_H

#pragma mark - mc_randn -

MC_TARGET_FUNC float mc_randnf(const float mu, const float stddev)
{
//!# Random number from Gaussian (normal) distribution with given mean and stddev.
	return mu + stddev * mc_randgf();
}

MC_TARGET_FUNC double mc_randnff(const float mu, const float stddev)
{
//!# Random number from Gaussian (normal) distribution with given mean and stddev.
	return mc_cast(double, mu) + mc_cast(double, stddev) * mc_randg();
}

MC_TARGET_FUNC double mc_randn(const double mu, const double stddev)
{
//!# Random number from Gaussian (normal) distribution with given mean and stddev.
	return mu + stddev * mc_randg();
}

MC_TARGET_FUNC long double mc_randnl(const long double mu, const long double stddev)
{
//!# Random number from Gaussian (normal) distribution with given mean and stddev.
	return mu + stddev * mc_randgl();
}

#endif /* !MC_RANDN_H */

/* EOF */