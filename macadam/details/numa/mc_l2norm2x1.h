//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_l2norm2x1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_hypot2.h>

#ifndef MC_NORM2X1_H
#define MC_NORM2X1_H

#pragma mark - mc_l2norm2x1 -

MC_TARGET_PROC float mc_l2norm2x1f(const int n, const int j, const float * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot2f(a[j], a[n + j]);
}

MC_TARGET_PROC double mc_l2norm2x1ff(const int n, const int j, const float * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot2(mc_cast(double, a[j]), mc_cast(double, a[n + j]));
}

MC_TARGET_PROC double mc_l2norm2x1(const int n, const int j, const double * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot2(a[j], a[n + j]);
}

MC_TARGET_PROC long double mc_l2norm2x1l(const int n, const int j, const long double * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot2l(a[j], a[n + j]);
}

#endif /* !MC_NORM2X1_H */

/* EOF */