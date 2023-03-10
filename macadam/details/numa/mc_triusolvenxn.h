//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_triusolvenxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zeros1xn.h>

#ifndef MC_TRIUSOLVENXN_H
#define MC_TRIUSOLVENXN_H

#pragma mark - MC_TRIUSOLVENXN_H -

MC_TARGET_FUNC int mc_triusolvenxnf(const int n, const float * u, const float * b, float * x)
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[n x n], b[n x 1], and x[1 x n].
	int i, j;
	float w;
	if (x != b) {
		for (i = 0; i < n; i++) {
			x[i] = b[i];
		}
	}
	for (j = (n - 1); j > 0; j--) {
		w = u[(n * j) + j];
		if (w != 0.0f) {
			x[j] = x[j] / w;
		}
		for (i = 0; i < j; i++) {
			x[i] = x[i] - x[j] * u[(n * i) + j];
		}
	}
	w = u[0];
	if (w != 0.0f) {
		x[0] = x[0] / w;
	}
	return 0;
}

MC_TARGET_FUNC int mc_triusolvenxnff(const int n, const float * u, const float * b, double * x)
{
//!# Solving the non-singular upper triangular system
//!# Ux=b, where u[n x n], b[n x 1], and x[1 x n].
	int i, j;
	double w;
	for (i = 0; i < n; i++) {
		x[i] = mc_cast(double, b[i]);
	}
	for (j = (n - 1); j > 0; j--) {
		w = mc_cast(double, u[(n * j) + j]);
		if (w != 0.0) {
			x[j] = x[j] / w;
		}
		for (i = 0; i < j; i++) {
			x[i] = x[i] - x[j] * mc_cast(double, u[(n * i) + j]);
		}
	}
	w = mc_cast(double, u[0]);
	if (w != 0.0) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_triusolvenxn(const int n, const double * u, const double * b, double * x)
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[n x n], b[n x 1], and x[1 x n].
	int i, j;
	double w;
	if (x != b) {
		for (i = 0; i < n; i++) {
			x[i] = b[i];
		}
	}
	for (j = (n - 1); j > 0; j--) {
		w = u[(n * j) + j];
		if (w != 0.0) {
			x[j] = x[j] / w;
		}
		for (i = 0; i < j; i++) {
			x[i] = x[i] - x[j] * u[(n * i) + j];
		}
	}
	w = u[0];
	if (w != 0.0) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_triusolvenxnl(const int n, const long double * u, const long double * b, long double * x)
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[n x n], b[n x 1], and x[1 x n].
	int i, j;
	long double w;
	if (x != b) {
		for (i = 0; i < n; i++) {
			x[i] = b[i];
		}
	}
	for (j = (n - 1); j > 0; j--) {
		w = u[(n * j) + j];
		if (w != 0.0L) {
			x[j] = x[j] / w;
		}
		for (i = 0; i < j; i++) {
			x[i] = x[i] - x[j] * u[(n * i) + j];
		}
	}
	w = u[0];
	if (w != 0.0L) {
		x[0] = x[0] / w;
	}

	return 0;
}

#endif /* !MC_TRIUSOLVENXN_H */

/* EOF */