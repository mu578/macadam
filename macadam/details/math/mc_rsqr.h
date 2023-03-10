//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rsqr.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>

#ifndef MC_RSQR_H
#define MC_RSQR_H

#pragma mark - mc_rsqr -

MC_TARGET_FUNC float mc_rsqrf(const float x)
{
	return 1.0f / mc_raise2f(x);
}

MC_TARGET_FUNC double mc_rsqr(const double x)
{
	return 1.0 / mc_raise2(x);
}

MC_TARGET_FUNC long double mc_rsqrl(const long double x)
{
	return 1.0L / mc_raise2l(x);
}

#endif /* !MC_RSQR_H */

/* EOF */