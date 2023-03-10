//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_xlogy.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log.h>

#ifndef MC_XLOGY_H
#define MC_XLOGY_H

#pragma mark - mc_xlogy -

MC_TARGET_FUNC float mc_xlogyf(const float x, const float y)
{
	return (x== 0 ? 0 : x * mc_logf(y));
}

MC_TARGET_FUNC double mc_xlogy(const double x, const double y)
{
	return (x== 0 ? 0 : x * mc_log(y));
}

MC_TARGET_FUNC long double mc_xlogyl(const long double x, const long double y)
{
	return (x== 0 ? 0 : x * mc_logl(y));
}

#endif /* !MC_XLOGY_H */

/* EOF */