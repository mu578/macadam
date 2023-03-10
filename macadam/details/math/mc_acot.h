//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_acot.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_atan.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ACOT_H
#define MC_ACOT_H

#pragma mark - mc_acot -

MC_TARGET_FUNC float mc_acotf(const float x)
{
	return mc_atanf(1.0f / x);
}

MC_TARGET_FUNC double mc_acot(const double x)
{
	return mc_atan(1.0 / x);
}

MC_TARGET_FUNC long double mc_acotl(const long double x)
{
	return mc_atanl(1.0L / x);
}

#endif /* !MC_ACOT_H */

/* EOF */