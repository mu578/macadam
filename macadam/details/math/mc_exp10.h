//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_exp10.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_pow.h>

#ifndef MC_EXP10_H
#define MC_EXP10_H

#pragma mark - mc_exp10 -

MC_TARGET_FUNC float mc_exp10f(const float x)
{
	return mc_powf(10.0f, x);
}

MC_TARGET_FUNC double mc_exp10(const double x)
{
	return mc_pow(10.0, x);
}

MC_TARGET_FUNC long double mc_exp10l(const long double x)
{
	return mc_powl(10.0L, x);
}

#endif /* !MC_EXP10_H */

/* EOF */