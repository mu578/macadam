//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_xlogp1y.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_xlog1py.h>

#ifndef MC_XLOGP1Y_H
#define MC_XLOGP1Y_H

#pragma mark - mc_xlogp1y -

MC_TARGET_FUNC float mc_xlogp1yf(const float x, const float y)
{
	return mc_xlog1pyf(x, y);
}

MC_TARGET_FUNC double mc_xlogp1y(const double x, const double y)
{
	return mc_xlog1py(x, y);
}

MC_TARGET_FUNC long double mc_xlogp1yl(const long double x, const long double y)
{
	return mc_xlog1pyl(x, y);
}

#endif /* !MC_XLOGP1Y_H */

/* EOF */