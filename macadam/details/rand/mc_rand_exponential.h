//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_exponential.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log1m.h>
#include <macadam/details/rand/mc_randu.h>

#ifndef MC_RAND_EXPONENTIAL_H
#define MC_RAND_EXPONENTIAL_H

#pragma mark - mc_rand_exponential -

MC_TARGET_FUNC float mc_rand_exponentialf(const float l)
{
//!# Exponential distribution generator.
	const float u = mc_randuf();
	const float s = 1.0f / (l > 0.0f ? l : 1.0f);
	return -mc_log1mf(u) * s;
}

MC_TARGET_FUNC double mc_rand_exponentialff(const float l)
{
//!# Exponential distribution generator.
	const double u = mc_randu();
	const double s = 1.0 / (mc_cast(double, l) > 0.0 ? mc_cast(double, l) : 1.0);
	return -mc_log1m(u) * s;
}

MC_TARGET_FUNC double mc_rand_exponential(const double l)
{
//!# Exponential distribution generator.
	const double u = mc_randu();
	const double s = 1.0 / (l > 0.0 ? l : 1.0);
	return -mc_log1m(u) * s;
}

MC_TARGET_FUNC long double mc_rand_exponentiall(const long double l)
{
//!# Exponential distribution generator.
	const long double u = mc_randul();
	const long double s = 1.0L / (l > 0.0L ? l : 1.0L);
	return -mc_log1ml(u) * s;
}

#endif /* !MC_RAND_EXPONENTIAL_H */

/* EOF */