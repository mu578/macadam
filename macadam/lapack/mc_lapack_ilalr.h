//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_ilalr.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_LAPACKE_ILALR_H
#define MC_LAPACKE_ILALR_H

#pragma mark - mc_lapack_ilaslr -

MC_TARGET_FUNC int mc_lapack_ilaslr(const int m, const int n, const float * a, const int lda)
{
	const float zero = 0.0f;

	int i, ilalr, j;

	if (m == 0) {
		mc_unused(lda);
		ilalr = m;
	} else if (mc_blas_matrix_at(a, lda, n, m, 1) != zero || mc_blas_matrix_at(a, lda, n, m, n) != zero) {
		ilalr = m;
	} else {
		ilalr = 0;
		for (j = 1; j <= n; ++j) {
			i = m;
			while (mc_blas_matrix_at(a, lda, n, mc_maxmag(i, 1), j) == zero && i >= 1) {
				i = i - 1;
			}
			ilalr = mc_maxmag(ilalr, i);
		}
	}
	return ilalr;
}

#pragma mark - mc_lapack_iladlr -

MC_TARGET_FUNC int mc_lapack_iladlr(const int m, const int n, const double * a, const int lda)
{
	const double zero = 0.0;

	int i, ilalr, j;

	if (m == 0) {
		mc_unused(lda);
		ilalr = m;
	} else if (mc_blas_matrix_at(a, lda, n, m, 1) != zero || mc_blas_matrix_at(a, lda, n, m, n) != zero) {
		ilalr = m;
	} else {
		ilalr = 0;
		for (j = 1; j <= n; ++j) {
			i = m;
			while (mc_blas_matrix_at(a, lda, n, mc_maxmag(i, 1), j) == zero && i >= 1) {
				i = i - 1;
			}
			ilalr = mc_maxmag(ilalr, i);
		}
	}
	return ilalr;
}

#pragma mark - mc_lapack_ilallr -

MC_TARGET_FUNC int mc_lapack_ilallr(const int m, const int n, const long double * a, const int lda)
{
	const long double zero = 0.0L;

	int i, ilalr, j;

	if (m == 0) {
		mc_unused(lda);
		ilalr = m;
	} else if (mc_blas_matrix_at(a, lda, n, m, 1) != zero || mc_blas_matrix_at(a, lda, n, m, n) != zero) {
		ilalr = m;
	} else {
		ilalr = 0;
		for (j = 1; j <= n; ++j) {
			i = m;
			while (mc_blas_matrix_at(a, lda, n, mc_maxmag(i, 1), j) == zero && i >= 1) {
				i = i - 1;
			}
			ilalr = mc_maxmag(ilalr, i);
		}
	}
	return ilalr;
}

#endif /* !MC_LAPACKE_ILALR_H */

/* EOF */