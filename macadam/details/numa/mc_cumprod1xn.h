//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cumprod1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_twoproduct.h>

#ifndef MC_CUMPROD1XN_H
#define MC_CUMPROD1XN_H

#pragma mark - mc_cumprod1xn -

MC_TARGET_FUNC int mc_cumprod1xnf(const int n, float * u, const float * x, const int f)
{
//!# Requires u[1 x n] and x[1 x n] where 1 < n. U and X may be
//!# the same. Returning the cumulative product of the elements.
//!# f=0: forward cumulative product.
//!# f=1: backward cumulative product.
	int i = 1;
	float a, b;
	if (n > 0) {
		u[(f == 1) ? 0 : (n - 1)] = x[(f == 1) ? 0 : (n - 1)];
		for (; i < n; i++) {
			const int k = (f == 1) ? (n - i - 1) : i;
			mc_twoproductf(u[k - 1], x[k], &a, &b);
			u[k]        = a + b;
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_cumprod1xnff(const int n, double * u, const float * x, const int f)
{
//!# Requires u[1 x n] and x[1 x n] where 1 < n. Returning
//!# the cumulative product of the elements.
//!# f=0: forward cumulative product.
//!# f=1: backward cumulative product.
	int i = 1;
	double a, b;
	if (n > 0) {
		u[(f == 1) ? 0 : (n - 1)] = x[(f == 1) ? 0 : (n - 1)];
		for (; i < n; i++) {
			const int k = (f == 1) ? (n - i - 1) : i;
			mc_twoproduct(u[k - 1], mc_cast(double, x[k]), &a, &b);
			u[k]        = a + b;
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_cumprod1xn(const int n, double * u, const double * x, const int f)
{
//!# Requires u[1 x n] and x[1 x n] where 1 < n. U and X may be
//!# the same. Returning the cumulative product of the elements.
//!# f=0: forward cumulative product.
//!# f=1: backward cumulative product.
	int i = 1;
	double a, b;
	if (n > 0) {
		u[(f == 1) ? 0 : (n - 1)] = x[(f == 1) ? 0 : (n - 1)];
		for (; i < n; i++) {
			const int k = (f == 1) ? (n - i - 1) : i;
			mc_twoproduct(u[k - 1], x[k], &a, &b);
			u[k]        = a + b;
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_cumprod1xnl(const int n, long double * u, const long double * x, const int f)
{
//!# Requires u[1 x n] and x[1 x n] where 1 < n. U and X may be
//!# the same. Returning the cumulative product of the elements.
//!# f=0: forward cumulative product.
//!# f=1: backward cumulative product.
	int i = 1;
	long double a, b;
	if (n > 0) {
		u[(f == 1) ? 0 : (n - 1)] = x[(f == 1) ? 0 : (n - 1)];
		for (; i < n; i++) {
			const int k = (f == 1) ? (n - i - 1) : i;
			mc_twoproductl(u[k - 1], x[k], &a, &b);
			u[k]        = a + b;
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_CUMPROD1XN_H */

/* EOF */