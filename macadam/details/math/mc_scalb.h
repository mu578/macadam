//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_scalb.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_scalbln.h>
#include <macadam/details/math/mc_nearbyint.h>

#ifndef MC_SCALB_H
#define MC_SCALB_H

#pragma mark - mc_scalb -

MC_TARGET_FUNC float mc_scalbf(const float x, const float y)
{
	return mc_scalblnf(x, mc_cast(const long, mc_nearbyintf(y)));
}

MC_TARGET_FUNC double mc_scalb(const double x, const double y)
{
	return mc_scalbln(x, mc_cast(const long, mc_nearbyint(y)));
}

MC_TARGET_FUNC long double mc_scalbl(const long double x, const long double y)
{
	return mc_scalblnl(x, mc_cast(const long, mc_nearbyintl(y)));
}

#endif /* !MC_SCALB_H */

/* EOF */