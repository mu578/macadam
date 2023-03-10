//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fix.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_ceil.h>
#include <macadam/details/math/mc_floor.h>

#ifndef MC_FIX_H
#define MC_FIX_H

#pragma mark - mc_fix -

MC_TARGET_FUNC float mc_fixf(const float x)
{
	return (x > 0.0f) ? mc_floorf(x) : ((x < 0.0f) ? mc_ceilf(x) : x);
}

MC_TARGET_FUNC double mc_fix(const double x)
{
	return (x > 0.0) ? mc_floor(x) : ((x < 0.0) ? mc_ceil(x) : x);
}

MC_TARGET_FUNC long double mc_fixl(const long double x)
{
	return (x > 0.0L) ? mc_floorl(x) : ((x < 0.0L) ? mc_ceill(x) : x);
}

#endif /* !MC_FIX_H */

/* EOF */