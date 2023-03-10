//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulatxmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/blas/mc_blas_gemv.h>

#ifndef MC_MULATXMXN_H
#define MC_MULATXMXN_H

#pragma mark - mc_mulatxmxn -

MC_TARGET_FUNC void mc_mulatxmxnf(const int m, const int n, float * MC_TARGET_RESTRICT b, const float * a, const float * x)
{
//!# Requires b[n x 1], a[m * n] and x[m x 1].
//!# b=a'*x
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm      = m;
	const int nn      = n;
	const int lda     = m;
	const int incx    = 1;
	const int incy    = 1;
	const float alpha = 1.0f;
	const float beta  = 0.0f;

	mc_blas_sgemv('T', mm, nn, alpha, a, lda, x, incx, beta, b, incy);
#	else
	int i, j = 0;
	float w;
	mc_zeros1xnf(n, b);
	for (; j < n; j++) {
		w = 0.0f;
		for (i = 0; i < m; i++) {
			w = w + (a[(n * i) + j] * x[i]);
		}
		b[j] = b[j] + w;
	}
#	endif
}

MC_TARGET_FUNC void mc_mulatxmxnff(const int m, const int n, double * MC_TARGET_RESTRICT b, const float * a, const float * x)
{
//!# Requires b[n x 1], a[m * n] and x[m x 1].
//!# b=a'*x
	int i, j = 0;
	double w;
	mc_zeros1xn(n, b);
	for (; j < n; j++) {
		w = 0.0;
		for (i = 0; i < m; i++) {
			w = w + (mc_cast(double, a[(n * i) + j]) * mc_cast(double, x[i]));
		}
		b[j] = b[j] + w;
	}
}

MC_TARGET_FUNC void mc_mulatxmxnfd(const int m, const int n, double * MC_TARGET_RESTRICT b, const float * a, const double * x)
{
//!# Requires b[n x 1], a[m * n] and x[m x 1].
//!# b=a'*x
	int i, j = 0;
	double w;
	mc_zeros1xn(n, b);
	for (; j < n; j++) {
		w = 0.0;
		for (i = 0; i < m; i++) {
			w = w + (mc_cast(double, a[(n * i) + j]) * x[i]);
		}
		b[j] = b[j] + w;
	}
}

MC_TARGET_FUNC void mc_mulatxmxndf(const int m, const int n, double * MC_TARGET_RESTRICT b, const double * a, const float * x)
{
//!# Requires b[n x 1], a[m * n] and x[m x 1].
//!# b=a'*x
	int i, j = 0;
	double w;
	mc_zeros1xn(n, b);
	for (; j < n; j++) {
		w = 0.0;
		for (i = 0; i < m; i++) {
			w = w + (a[(n * i) + j] * mc_cast(double, x[i]));
		}
		b[j] = b[j] + w;
	}
}

MC_TARGET_FUNC void mc_mulatxmxn(const int m, const int n, double * MC_TARGET_RESTRICT b, const double * a, const double * x)
{
//!# Requires b[n x 1], a[m * n] and x[m x 1].
//!# b=a'*x
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm       = m;
	const int nn       = n;
	const int lda      = m;
	const int incx     = 1;
	const int incy     = 1;
	const double alpha = 1.0;
	const double beta  = 0.0;

	mc_blas_dgemv('T', mm, nn, alpha, a, lda, x, incx, beta, b, incy);
#	else
	int i, j = 0;
	double w;
	mc_zeros1xn(n, b);
	for (; j < n; j++) {
		w = 0.0;
		for (i = 0; i < m; i++) {
			w = w + (a[(n * i) + j] * x[i]);
		}
		b[j] = b[j] + w;
	}
#	endif
}

MC_TARGET_FUNC void mc_mulatxmxnl(const int m, const int n, long double * MC_TARGET_RESTRICT b, const long double * a, const long double * x)
{
//!# Requires b[n x 1], a[m * n] and x[m x 1].
//!# b=a'*x
#	if !MC_TARGET_EMBEDDED && MC_TARGET_BLAS_USE_CLAYOUT

	const int mm            = m;
	const int nn            = n;
	const int lda           = m;
	const int incx          = 1;
	const int incy          = 1;
	const long double alpha = 1.0L;
	const long double beta  = 0.0L;

	mc_blas_lgemv('T', mm, nn, alpha, a, lda, x, incx, beta, b, incy);
#	else
	int i, j = 0;
	long double w;
	mc_zeros1xnl(n, b);
	for (; j < n; j++) {
		w = 0.0L;
		for (i = 0; i < m; i++) {
			w = w + (a[(n * i) + j] * x[i]);
		}
		b[j] = b[j] + w;
	}
#	endif
}

#endif /* !MC_MULATXMXN_H */

/* EOF */