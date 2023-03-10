//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_diff1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_DIFF1XN_H
#define MC_DIFF1XN_H

#pragma mark - mc_diff1xn -

MC_TARGET_FUNC int mc_diff1xnf(const int n, float * y, const float * x)
{
//!# Requires y[n - 1] and x[n] where 1 < n. Y and X may be the same.
//!# Calculating differences between adjacent elements of vector `x`.
//!# Returns the number of new elements in y vector.
	int i = 1;
	for (; i < n; i++) {
		y[i - 1] = x[i] - x[i - 1];
	}
	return n - 1;
}

MC_TARGET_FUNC int mc_diff1xnff(const int n, double * y, const float * x)
{
//!# Requires y[n - 1] and x[n] where 1 < n. Calculating differences
//!# between adjacent elements of vector `x`. Returns the number of
//!# elements in y vectors.
	int i = 1;
	for (; i < n; i++) {
		y[i - 1] = mc_cast(double, x[i]) - mc_cast(double, x[i - 1]);
	}
	return n - 1;
}

MC_TARGET_FUNC int mc_diff1xn(const int n, double * y, const double * x)
{
//!# Requires y[n - 1] and x[n] where 1 < n. Y and X may be the same.
//!# Calculating differences between adjacent elements of vector `x`.
//!# Returns the number of new elements in y vector.
	int i = 1;
	for (; i < n; i++) {
		y[i - 1] = x[i] - x[i - 1];
	}
	return n - 1;
}

MC_TARGET_FUNC int mc_diff1xnl(const int n, long double * y, const long double * x)
{
//!# Requires y[n - 1] and x[n] where 1 < n. Y and X may be the same.
//!# Calculating differences between adjacent elements of vector `x`.
//!# Returns the number of new elements in y vector.
	int i = 1;
	for (; i < n; i++) {
		y[i - 1] = x[i] - x[i - 1];
	}
	return n - 1;
}

#endif /* !MC_DIFF1XN_H */

/* EOF */