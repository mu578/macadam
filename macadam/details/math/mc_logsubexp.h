//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logsubexp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_logdiffexp.h>

#ifndef MC_LOGSUBEXP_H
#define MC_LOGSUBEXP_H

#pragma mark - mc_logsubexp -

MC_TARGET_FUNC float mc_logsubexpf(float x, float y)
{
	if (x == y) {
		return MCK_INFN;
	}
	if (x < y) {
		const float w = x;
		x             = y;
		y             = w;
	}
	if (!mc_isinf(y) && y > MCK_INFN) {
		return mc_logdiffexpf(x, y);
	}
	return x;
}

MC_TARGET_FUNC double mc_logsubexp(double x, double y)
{
	if (x == y) {
		return MCK_INFN;
	}
	if (x < y) {
		const double w = x;
		x              = y;
		y              = w;
	}
	if (!mc_isinf(y) && y > MCK_INFN) {
		return mc_logdiffexp(x, y);
	}
	return x;
}

MC_TARGET_FUNC long double mc_logsubexpl(long double x, long double y)
{
	if (x == y) {
		return MCK_INFN;
	}
	if (x < y) {
		const long double w = x;
		x                   = y;
		y                   = w;
	}
	if (!mc_isinf(y) && y > MCK_INFN) {
		return mc_logdiffexpl(x, y);
	}
	return x;
}

#endif /* !MC_LOGSUBEXP_H */

/* EOF */