//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_linspace1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_LINSPACE1XN
#define MC_LINSPACE1XN

#pragma mark - mc_linspace1xn -

MC_TARGET_FUNC int mc_linspace1xnf(const int n, float * x, float x1, float x2)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a linearly spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	float step;
	if (n > 0) {
		if (n < 2) {
			x[0] = x2;
		} else if (n < 3) {
			x[0] = x1;
			x[1] = x2;
		} else {
			step       = (x2 - x1) / mc_cast(const float, (n - 1));
			x[0]       = x1;
			x[(n - 1)] = x2;
			for (; i < (n - 1); i++) {
				x[i] = x1 + (mc_cast(float, i) * step);
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_linspace1xnff(const int n, double * x, float x1, float x2)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a linearly spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	double step, x1d, x2d;
	if (n > 0) {
		x1d = mc_cast(double, x1);
		x2d = mc_cast(double, x2);
		if (n < 2) {
			x[0] = x2d;
		} else if (n < 3) {
			x[0] = x1d;
			x[1] = x2d;
		} else {
			step       = (x2d - x1d) / mc_cast(const double, (n - 1));
			x[0]       = x1d;
			x[(n - 1)] = x2d;
			for (; i < (n - 1); i++) {
				x[i] = x1d +  (mc_cast(double, i) * step);
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_linspace1xn(const int n, double * x, double x1, double x2)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a linearly spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	double step;
	if (n > 0) {
		if (n < 2) {
			x[0] = x2;
		} else if (n < 3) {
			x[0] = x1;
			x[1] = x2;
		} else {
			step       = (x2 - x1) / mc_cast(const double, (n - 1));
			x[0]       = x1;
			x[(n - 1)] = x2;
			for (; i < (n - 1); i++) {
				x[i] = x1 + (mc_cast(double, i) * step);
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_linspace1xnl(const int n, long double * x, long double x1, long double x2)
{
//!# Requires x[n] where 1 < n. Draws a linspace: generates a linearly spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	long double step;
	if (n > 0) {
		if (n < 2) {
			x[0] = x2;
		} else if (n < 3) {
			x[0] = x1;
			x[1] = x2;
		} else {
			step       = (x2 - x1) / mc_cast(const long double, (n - 1));
			x[0]       = x1;
			x[(n - 1)] = x2;
			for (; i < (n - 1); i++) {
				x[i] = x1 + (mc_cast(long double, i) * step);
			}
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_LINSPACE1XN */

/* EOF */