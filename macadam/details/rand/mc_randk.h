//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_randk.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/rand/mc_randi.h>

#ifndef MC_RANDK_H
#define MC_RANDK_H

#pragma mark - mc_randk -

MC_TARGET_PROC float mc_randkf(void)
{
//!# Uniform distribution range [0, 1]. Using an arbitrary
//!# pair of integers for generating sample over field k.
	const unsigned int k = 0x1000001U;
	const unsigned int a = mc_randi() % k;
	const unsigned int b = mc_randi() % k;
	return mc_cast_expr(float, a & 1 ? b : (b & 1 ? a : b)) / mc_cast(const float, k);
}

MC_TARGET_PROC double mc_randk(void)
{
//!# Uniform distribution range [0, 1]. Using an arbitrary
//!# pair of integers for generating sample over field k.
	const unsigned int k = MCLIMITS_RANDMAX;
	const unsigned int a = mc_randi() % k;
	const unsigned int b = mc_randi() % k;
	return mc_cast_expr(double, a & 1 ? b : (b & 1 ? a : b)) / mc_cast(const double, k);
}

MC_TARGET_PROC long double mc_randkl(void)
{
//!# Uniform distribution range [0, 1]. Using an arbitrary
//!# pair of integers for generating sample over field k.
	const unsigned int k = MCLIMITS_RANDMAX;
	const unsigned int a = mc_randi() % k;
	const unsigned int b = mc_randi() % k;
	return mc_cast_expr(long double, a & 1 ? b : (b & 1 ? a : b)) / mc_cast(const long double, k);
}

#endif /* !MC_RANDK_H */

/* EOF */