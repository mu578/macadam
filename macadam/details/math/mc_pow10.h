//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_pow10.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp10.h>

#ifndef MC_POW10_H
#define MC_POW10_H

#pragma mark - mc_pow10 -

MC_TARGET_FUNC float mc_pow10f(const float x)
{
	return mc_exp10f(x);
}

MC_TARGET_FUNC double mc_pow10(const double x)
{
	return mc_exp10(x);
}

MC_TARGET_FUNC long double mc_pow10l(const long double x)
{
	return mc_exp10l(x);
}

#endif /* !MC_POW10_H */

/* EOF */