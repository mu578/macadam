//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_weibull.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/rand/mc_rand_exponential.h>

#ifndef MC_RAND_WEIBULL_H
#define MC_RAND_WEIBULL_H

#pragma mark - mc_rand_weibull -

MC_TARGET_FUNC float mc_rand_weibullf(const float a)
{
//!# Weibull distribution generator. a=alpha=shape.
	return a != 0.0f ? (a == 1.0f ? mc_rand_exponentialf(1.0f) : mc_powf(mc_rand_exponentialf(1.0f), 1.0f / a)) : 0.0f;
}

MC_TARGET_FUNC double mc_rand_weibullff(const float a)
{
//!# Weibull distribution generator. a=alpha=shape.
	return a != 0.0f ? (a == 1.0f ? mc_rand_exponential(1.0) : mc_pow(mc_rand_exponential(1.0), 1.0 / mc_cast(const double, a))) : 0.0;
}

MC_TARGET_FUNC double mc_rand_weibull(const double a)
{
//!# Weibull distribution generator. a=alpha=shape.
	return a != 0.0 ? (a == 1.0 ? mc_rand_exponential(1.0) : mc_pow(mc_rand_exponential(1.0), 1.0 / a)) : 0.0;
}

MC_TARGET_FUNC long double mc_rand_weibulll(const long double a)
{
//!# Weibull distribution generator. a=alpha=shape.
	return a != 0.0L ? (a == 1.0L ? mc_rand_exponentiall(1.0L) : mc_powl(mc_rand_exponentiall(1.0), 1.0L / a)) : 0.0L;
}

#endif /* !MC_RAND_WEIBULL_H */

/* EOF */