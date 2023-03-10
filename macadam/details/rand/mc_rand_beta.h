//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_beta.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/rand/mc_rand_gamma.h>

#ifndef MC_RAND_BETA_H
#define MC_RAND_BETA_H

#pragma mark - mc_rand_beta -

MC_TARGET_FUNC float mc_rand_betaf(const float a, const float b)
{
//!# Beta distribution generator.
	const float x = mc_rand_gammaf(a, 1.0f);
	const float y = mc_rand_gammaf(b, 1.0f);
	const float z = x + y;
	return z != 0.0f ? x / z : 0.0f;
}

MC_TARGET_FUNC double mc_rand_betaff(const float a, const float b)
{
//!# Beta distribution generator.
	const double x = mc_rand_gammaff(a, 1.0f);
	const double y = mc_rand_gammaff(b, 1.0f);
	const double z = x + y;
	return z != 0.0 ? x / z : 0.0;
}

MC_TARGET_FUNC double mc_rand_beta(const double a, const double b)
{
//!# Beta distribution generator.
	const double x = mc_rand_gamma(a, 1.0);
	const double y = mc_rand_gamma(b, 1.0);
	const double z = x + y;
	return z != 0.0 ? x / z : 0.0;
}

MC_TARGET_FUNC long double mc_rand_betal(const long double a, const long double b)
{
//!# Beta distribution generator.
	const long double x = mc_rand_gammal(a, 1.0L);
	const long double y = mc_rand_gammal(b, 1.0L);
	const long double z = x + y;
	return z != 0.0L ? x / z : 0.0L;
}

#endif /* !MC_RAND_BETA_H */

/* EOF */