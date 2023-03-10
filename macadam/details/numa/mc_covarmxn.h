//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_covarmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_meanmx1.h>

#ifndef MC_COVARMXN_H
#define MC_COVARMXN_H

#pragma mark - mc_covarmxn -

MC_TARGET_FUNC int mc_covarmxnf(const int m, const int n, float * MC_TARGET_RESTRICT c, const float * a, const int b)
{
//!# Requires c[n x n] and a[m x n] where 1 < n <= m. Estimating a covariance matrix C
//!# from A observations. A is a matrix whose columns represent random variables and
//!# whose rows represent observations. C is the covariance matrix with the corresponding
//!# column variances along the diagonal.
	int i, j = 0, k;
	float mj, mk;
	if (1 < n && n <= m) {
		for (; j < n; j++) {
			mj = mc_meanmx1f(m, n, j, a, 0, 5);
			for (k = 0; k < n; k++) {
				mk             = mc_meanmx1f(m, n, k, a, 0, 5);
				c[(n * j) + k] = 0.0f;
				for (i = 0; i < m; i++) {
					c[(n * j) + k] = c[(n * j) + k] + ((a[(n * i) + j] - mj) * (a[(n * i) + k] - mk));
				}
				c[(n * j) + k] = c[(n * j) + k] / mc_cast(const float, (b ? m - 1 : m));
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_covarmxnff(const int m, const int n, double * c, const float * a, const int b)
{
//!# Requires c[n x n] and a[m x n] where 1 < n <= m. Estimating a covariance matrix C
//!# from A observations. A is a matrix whose columns represent random variables and
//!# whose rows represent observations. C is the covariance matrix with the corresponding
//!# column variances along the diagonal.
	int i, j = 0, k;
	double mj, mk;
	if (1 < n && n <= m) {
		for (; j < n; j++) {
			mj = mc_cast(double, mc_meanmx1f(m, n, j, a, 0, 5));
			for (k = 0; k < n; k++) {
				mk             = mc_cast(double, mc_meanmx1f(m, n, k, a, 0, 5));
				c[(n * j) + k] = 0.0;
				for (i = 0; i < m; i++) {
					c[(n * j) + k] = c[(n * j) + k] + ((mc_cast(double, a[(n * i) + j]) - mj) * (mc_cast(double, a[(n * i) + k]) - mk));
				}
				c[(n * j) + k] = c[(n * j) + k] / mc_cast(const double, (b ? m - 1 : m));
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_covarmxn(const int m, const int n, double * MC_TARGET_RESTRICT c, const double * a, const int b)
{
//!# Requires c[n x n] and a[m x n] where 1 < n <= m. Estimating a covariance matrix C
//!# from A observations. A is a matrix whose columns represent random variables and
//!# whose rows represent observations. C is the covariance matrix with the corresponding
//!# column variances along the diagonal.
	int i, j = 0, k;
	double mj, mk;
	if (1 < n && n <= m) {
		for (; j < n; j++) {
			mj = mc_meanmx1(m, n, j, a, 0, 5);
			for (k = 0; k < n; k++) {
				mk             = mc_meanmx1(m, n, k, a, 0, 5);
				c[(n * j) + k] = 0.0;
				for (i = 0; i < m; i++) {
					c[(n * j) + k] = c[(n * j) + k] + ((a[(n * i) + j] - mj) * (a[(n * i) + k] - mk));
				}
				c[(n * j) + k] = c[(n * j) + k] / mc_cast(const double, (b ? m - 1 : m));
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_covarmxnl(const int m, const int n, long double * MC_TARGET_RESTRICT c, const long double * a, const int b)
{
//!# Requires c[n x n] and a[m x n] where 1 < n <= m. Estimating a covariance matrix C
//!# from A observations. A is a matrix whose columns represent random variables and
//!# whose rows represent observations. C is the covariance matrix with the corresponding
//!# column variances along the diagonal.
	int i, j = 0, k;
	long double mj, mk;
	if (1 < n && n <= m) {
		for (; j < n; j++) {
			mj = mc_meanmx1l(m, n, j, a, 0, 5);
			for (k = 0; k < n; k++) {
				mk             = mc_meanmx1l(m, n, k, a, 0, 5);
				c[(n * j) + k] = 0.0L;
				for (i = 0; i < m; i++) {
					c[(n * j) + k] = c[(n * j) + k] + ((a[(n * i) + j] - mj) * (a[(n * i) + k] - mk));
				}
				c[(n * j) + k] = c[(n * j) + k] / mc_cast(const long double, (b ? m - 1 : m));
			}
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_COVARMXN_H */

/* EOF */