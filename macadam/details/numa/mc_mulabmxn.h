//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulabmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>

#ifndef MC_MULABMXN_H
#define MC_MULABMXN_H

#pragma mark - mc_mulabmxn -

MC_TARGET_FUNC void mc_mulabmxnf(const int m, const int n, const int p, float * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# c=a*b. Producing c[m x p]=a[m x n] * b[n x p].

#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm      = m;
	const int nn      = p;
	const int kk      = n;
	const int lda     = m;
	const int ldb     = n;
	const int ldc     = m;
	const float alpha = 1.0f;
	const float beta  = 0.0f;

	mc_sgemm('N', 'N', mm, nn, kk, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	int i = 0, j, k;
	for (; i < m; i++) {
		for (j = 0; j < p; j++) {
			c[(p * i) + j] = 0.0f;
			for (k = 0; k < n; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * i) + k] * b[(p * k) + j]);
			}
		}
	}
#	endif
}

MC_TARGET_FUNC void mc_mulabmxnff(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# c=a*b. Producing c[m x p]=a[m x n] * b[n x p].
	int i = 0, j, k;
	for (; i < m; i++) {
		for (j = 0; j < p; j++) {
			c[(p * i) + j] = 0.0f;
			for (k = 0; k < n; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (mc_cast(double, a[(n * i) + k]) * mc_cast(double, b[(p * k) + j]));
			}
		}
	}
}

MC_TARGET_FUNC void mc_mulabmxnfd(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const float * a, const double * b)
{
//!# c=a*b. Producing c[m x p]=a[m x n] * b[n x p].
	int i = 0, j, k;
	for (; i < m; i++) {
		for (j = 0; j < p; j++) {
			c[(p * i) + j] = 0.0f;
			for (k = 0; k < n; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (mc_cast(double, a[(n * i) + k]) * b[(p * k) + j]);
			}
		}
	}
}

MC_TARGET_FUNC void mc_mulabmxndf(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const double * a, const float * b)
{
//!# c=a*b. Producing c[m x p]=a[m x n] * b[n x p].
	int i = 0, j, k;
	for (; i < m; i++) {
		for (j = 0; j < p; j++) {
			c[(p * i) + j] = 0.0f;
			for (k = 0; k < n; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * i) + k] * mc_cast(double, b[(p * k) + j]));
			}
		}
	}
}

MC_TARGET_FUNC void mc_mulabmxn(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const double * a, const double * b)
{
//!# c=a*b. Producing c[m x p]=a[m x n] * b[n x p].
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm       = m;
	const int nn       = p;
	const int kk       = n;
	const int lda      = m;
	const int ldb      = n;
	const int ldc      = m;
	const double alpha = 1.0;
	const double beta  = 0.0;

	mc_dgemm('N', 'N', mm, nn, kk, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	int i = 0, j, k;
	for (; i < m; i++) {
		for (j = 0; j < p; j++) {
			c[(p * i) + j] = 0.0;
			for (k = 0; k < n; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * i) + k] * b[(p * k) + j]);
			}
		}
	}
#	endif
}

MC_TARGET_FUNC void mc_mulabmxnl(const int m, const int n, const int p, long double * MC_TARGET_RESTRICT c, const long double * a, const long double * b)
{
//!# c=a*b. Producing c[m x p]=a[m x n] * b[n x p].
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm            = m;
	const int nn            = p;
	const int kk            = n;
	const int lda           = m;
	const int ldb           = n;
	const int ldc           = m;
	const long double alpha = 1.0L;
	const long double beta  = 0.0L;

	mc_lgemm('N', 'N', mm, nn, kk, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	int i = 0, j, k;
	for (; i < m; i++) {
		for (j = 0; j < p; j++) {
			c[(p * i) + j] = 0.0L;
			for (k = 0; k < n; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * i) + k] * b[(p * k) + j]);
			}
		}
	}
#	endif
}

#endif /* !MC_MULABMXN_H */

/* EOF */