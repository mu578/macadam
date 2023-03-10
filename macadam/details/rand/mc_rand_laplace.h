//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_laplace.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/rand/mc_rand_exponential.h>
#include <macadam/details/rand/mc_randu.h>

#ifndef MC_RAND_LAPLACE_H
#define MC_RAND_LAPLACE_H

#pragma mark - mc_rand_laplace -

MC_TARGET_FUNC float mc_rand_laplacef(const float l)
{
//!# Laplace distribution generator.
	const float u = mc_randuf() < 0.5f ? -1.0f : 1.0f;
	const float e = mc_rand_exponentialf(MCK_KF(MCK_1_SQRT2) * l);
	return u * e;
}

MC_TARGET_FUNC double mc_rand_laplaceff(const float l)
{
//!# Laplace distribution generator.
	const double u = mc_randu() < 0.5 ? -1.0 : 1.0;
	const double e = mc_rand_exponential(MCK_K(MCK_1_SQRT2) * mc_cast(double, l));
	return u * e;
}

MC_TARGET_FUNC double mc_rand_laplace(const double l)
{
//!# Laplace distribution generator.
	const double u = mc_randu() < 0.5 ? -1.0 : 1.0;
	const double e = mc_rand_exponential(MCK_K(MCK_1_SQRT2) * l);
	return u * e;
}

MC_TARGET_FUNC long double mc_rand_laplacel(const long double l)
{
//!# Laplace distribution generator.
	const long double u = mc_randul() < 0.5L ? -1.0L : 1.0L;
	const long double e = mc_rand_exponentiall(MCK_KL(MCK_1_SQRT2) * l);
	return u * e;
}

#endif /* !MC_RAND_LAPLACE_H */

/* EOF */