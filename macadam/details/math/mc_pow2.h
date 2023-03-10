//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_pow2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_floor.h>
#include <macadam/details/math/mc_exp2.h>
#include <macadam/details/math/mc_exp2i.h>

#ifndef MC_POW2_H
#define MC_POW2_H

#pragma mark - mc_pow2 -

MC_TARGET_FUNC float mc_pow2f(const float x)
{
	if (mc_floorf(x) == x) {
		return mc_exp2if(mc_cast(const int, x));
	}
	return mc_exp2f(x);
}

MC_TARGET_FUNC double mc_pow2(const double x)
{
	if (mc_floor(x) == x) {
		return mc_exp2i(mc_cast(const int, x));
	}
	return mc_exp2(x);
}

MC_TARGET_FUNC long double mc_pow2l(const long double x)
{
	if (mc_floorl(x) == x) {
		return mc_exp2il(mc_cast(const int, x));
	}
	return mc_exp2l(x);
}

#endif /* !MC_POW2_H */

/* EOF */