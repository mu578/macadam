//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_arange1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ARANGE1XN
#define MC_ARANGE1XN

#pragma mark - mc_arange1xn -

MC_TARGET_FUNC int mc_arangef(const int n, float * x, float x1, float x2, float xs)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a evenly
//!# spaced values within a given interval, values are generated within
//!# the half-open interval [x1, x2[. Returning usable length of vector X.
	int i   = 0;
	float w = x1;
	if (n > 0) {
		for (; (xs >= 0.0f) ? w < x2 : w > x2 && i < n; w += xs, i++) {
			x[i] = w;
		}
		return i;
	}
	return -1;
}

MC_TARGET_FUNC int mc_arangeff(const int n, double * x, float x1, float x2, float xs)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a evenly
//!# spaced values within a given interval, values are generated within
//!# the half-open interval [x1, x2[. Returning usable length of vector X.
	int i    = 0;
	double w = mc_cast(double, x1);
	if (n > 0) {
		for (; (xs >= 0.0f) ? w < mc_cast(double, x2) : w > mc_cast(double, x2) && i < n; w += mc_cast(double, xs), i++) {
			x[i] = w;
		}
		return i;
	}
	return -1;
}

MC_TARGET_FUNC int mc_arange(const int n, double * x, double x1, double x2, double xs)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a evenly
//!# spaced values within a given interval, values are generated within
//!# the half-open interval [x1, x2[. Returning usable length of vector X.
	int i    = 0;
	double w = x1;
	if (n > 0) {
		for (; (xs >= 0.0) ? w < x2 : w > x2 && i < n; w += xs, i++) {
			x[i] = w;
		}
		return i;
	}
	return -1;
}

MC_TARGET_FUNC int mc_arangel(const int n, long double * x, long double x1, long double x2, long double xs)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a evenly
//!# spaced values within a given interval, values are generated within
//!# the half-open interval [x1, x2[. Returning usable length of vector X.
	int i         = 0;
	long double w = x1;
	if (n > 0) {
		for (; (xs >= 0.0L) ? w < x2 : w > x2 && i < n; w += xs, i++) {
			x[i] = w;
		}
		return i;
	}
	return -1;
}

#endif /* !MC_ARANGE1XN */

/* EOF */