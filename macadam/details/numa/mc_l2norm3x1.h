//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_l2norm3x1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_hypot3.h>

#ifndef MC_NORM3X1_H
#define MC_NORM3X1_H

#pragma mark - mc_l2norm3x1 -

MC_TARGET_PROC float mc_l2norm3x1f(const int n, const int j, const float * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot3f(a[j], a[n + j], a[(n * 2) + j]);
}

MC_TARGET_PROC double mc_l2norm3x1ff(const int n, const int j, const float * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot3(mc_cast(double, a[j]), mc_cast(double, a[n + j]), mc_cast(double, a[(n * 2) + j]));
}

MC_TARGET_PROC double mc_l2norm3x1(const int n, const int j, const double * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot3(a[j], a[n + j], a[(n * 2) + j]);
}

MC_TARGET_PROC long double mc_l2norm3x1l(const int n, const int j, const long double * a)
{
	if (n < 1) {
		return MCK_NAN;
	}
	return mc_hypot3l(a[j], a[n + j], a[(n * 2) + j]);
}

#endif /* !MC_NORM3X1_H */

/* EOF */