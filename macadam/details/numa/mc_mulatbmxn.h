//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulatbmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>

#ifndef MC_MULATBMXN_H
#define MC_MULATBMXN_H

#pragma mark - mc_mulatbmxn -

MC_TARGET_FUNC void mc_mulatbmxnf(const int m, const int n, const int p, float * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# c=a'*b. Producing c[m x p]=a[m x n] * b[n x p].
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm      = n;
	const int nn      = p;
	const int kk      = m;
	const int lda     = m;
	const int ldb     = m;
	const int ldc     = n;
	const float alpha = 1.0f;
	const float beta  = 0.0f;

	mc_sgemm('T', 'N', mm, nn, kk, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	int i, j = 0, k;
	for (; j < p; j++) {
		for (i = 0; i < n; i++) {
			c[(p * i) + j] = 0.0f;
			for (k = 0; k < m; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * k) + i] * b[(p * k) + j]);
			}
		}
	}
#	endif
}

MC_TARGET_FUNC void mc_mulatbmxnff(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# c=a'*b. Producing c[m x p]=a[m x n] * b[n x p].
	int i, j = 0, k;
	for (; j < p; j++) {
		for (i = 0; i < n; i++) {
			c[(p * i) + j] = 0.0;
			for (k = 0; k < m; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (mc_cast(double, a[(n * k) + i]) * mc_cast(double, b[(p * k) + j]));
			}
		}
	}
}

MC_TARGET_FUNC void mc_mulatbmxnfd(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const float * a, const double * b)
{
//!# c=a'*b. Producing c[m x p]=a[m x n] * b[n x p].
	int i, j = 0, k;
	for (; j < p; j++) {
		for (i = 0; i < n; i++) {
			c[(p * i) + j] = 0.0;
			for (k = 0; k < m; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (mc_cast(double, a[(n * k) + i]) * b[(p * k) + j]);
			}
		}
	}
}

MC_TARGET_FUNC void mc_mulatbmxndf(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const double * a, const float * b)
{
//!# c=a'*b. Producing c[m x p]=a[m x n] * b[n x p].
	int i, j = 0, k;
	for (; j < p; j++) {
		for (i = 0; i < n; i++) {
			c[(p * i) + j] = 0.0;
			for (k = 0; k < m; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * k) + i] * mc_cast(double, b[(p * k) + j]));
			}
		}
	}
}

MC_TARGET_FUNC void mc_mulatbmxn(const int m, const int n, const int p, double * MC_TARGET_RESTRICT c, const double * a, const double * b)
{
//!# c=a'*b. Producing c[m x p]=a[m x n] * b[n x p].
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm       = n;
	const int nn       = p;
	const int kk       = m;
	const int lda      = m;
	const int ldb      = m;
	const int ldc      = n;
	const double alpha = 1.0;
	const double beta  = 0.0;

	mc_dgemm('T', 'N', mm, nn, kk, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	int i, j = 0, k;
	for (; j < p; j++) {
		for (i = 0; i < n; i++) {
			c[(p * i) + j] = 0.0;
			for (k = 0; k < m; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * k) + i] * b[(p * k) + j]);
			}
		}
	}
#	endif
}

MC_TARGET_FUNC void mc_mulatbmxnl(const int m, const int n, const int p, long double * MC_TARGET_RESTRICT c, const long double * a, const long double * b)
{
//!# c=a'*b. Producing c[m x p]=a[m x n] * b[n x p].
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm            = n;
	const int nn            = p;
	const int kk            = m;
	const int lda           = m;
	const int ldb           = m;
	const int ldc           = n;
	const long double alpha = 1.0L;
	const long double beta  = 0.0L;

	mc_lgemm('T', 'N', mm, nn, kk, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	int i, j = 0, k;
	for (; j < p; j++) {
		for (i = 0; i < n; i++) {
			c[(p * i) + j] = 0.0L;
			for (k = 0; k < m; k++) {
				c[(p * i) + j] = c[(p * i) + j] + (a[(n * k) + i] * b[(p * k) + j]);
			}
		}
	}
#	endif
}

#endif /* !MC_MULATBMXN_H */

/* EOF */