//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trilsolvenxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zeros1xn.h>

#ifndef MC_TRILSOLVENXN_H
#define MC_TRILSOLVENXN_H

#pragma mark - MC_TRILSOLVENXN_H -

MC_TARGET_FUNC int mc_trilsolvenxnf(const int n, const float * l, const float * b, float * x)
{
//!# B and X may be the same. Solving the non-singular lower
//!# triangular system Lx=b, where l[n x n], b[n x 1], and x[1 x n].
	int i, j;
	float w;
	if (x != b) {
		for (i = 0; i < n; i++) {
			x[i] = b[i];
		}
	}
	for (j = 0; j < (n - 1); j++) {
		w    = l[(n * j) + j];
		if (w != 0.0f) {
			x[j] = x[j] / w;
		}
		for (i = j + 1; i < n; i++) {
			x[i] = x[i] - l[(n * i) + j] * x[j];
		}
	}
	w = l[(n * n) - 1];
	if (w != 0.0f) {
		x[n - 1] = x[n - 1] / w;
	}
	return 0;
}

MC_TARGET_FUNC int mc_trilsolvenxnff(const int n, const float * l, const float * b, double * x)
{
//!# Solving the non-singular lower triangular system
//!# Lx=b, where l[n x n], b[n x 1], and x[1 x n].
	int i, j;
	double w;
	for (i = 0; i < n; i++) {
		x[i] = mc_cast(double, b[i]);
	}
	for (j = 0; j < (n - 1); j++) {
		w    = mc_cast(double, l[(n * j) + j]);
		if (w != 0.0) {
			x[j] = x[j] / w;
		}
		for (i = j + 1; i < n; i++) {
			x[i] = x[i] - l[(n * i) + j] * x[j];
		}
	}
	w = mc_cast(double, l[(n * n) - 1]);
	if (w != 0.0) {
		x[n - 1] = x[n - 1] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_trilsolvenxn(const int n, const double * l, const double * b, double * x)
{
//!# B and X may be the same. Solving the non-singular lower
//!# triangular system Lx=b, where l[n x n], b[n x 1], and x[1 x n].
	int i, j;
	double w;
	if (x != b) {
		for (i = 0; i < n; i++) {
			x[i] = b[i];
		}
	}
	for (j = 0; j < (n - 1); j++) {
		w    = l[(n * j) + j];
		if (w != 0.0) {
			x[j] = x[j] / w;
		}
		for (i = j + 1; i < n; i++) {
			x[i] = x[i] - l[(n * i) + j] * x[j];
		}
	}
	w = l[(n * n) - 1];
	if (w != 0.0) {
		x[n - 1] = x[n - 1] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_trilsolvenxnl(const int n, const long double * l, const long double * b, long double * x)
{
//!# B and X may be the same. Solving the non-singular lower
//!# triangular system Lx=b, where l[n x n], b[n x 1], and x[1 x n].
	int i, j;
	long double w;
	if (x != b) {
		for (i = 0; i < n; i++) {
			x[i] = b[i];
		}
	}
	for (j = 0; j < (n - 1); j++) {
		w    = l[(n * j) + j];
		if (w != 0.0L) {
			x[j] = x[j] / w;
		}
		for (i = j + 1; i < n; i++) {
			x[i] = x[i] - l[(n * i) + j] * x[j];
		}
	}
	w = l[(n * n) - 1];
	if (w != 0.0L) {
		x[n - 1] = x[n - 1] / w;
	}

	return 0;
}

#endif /* !MC_TRILSOLVENXN_H */

/* EOF */