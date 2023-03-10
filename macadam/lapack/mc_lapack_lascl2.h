//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lascl2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/lapack/mc_lapack_lartgp.h>
#include <macadam/details/math/mc_fabs.h>

#ifndef MC_LAPACKE_LASCL2_H
#define MC_LAPACKE_LASCL2_H

#pragma mark - mc_lapack_slascl2 -

MC_TARGET_FUNC void mc_lapack_slascl2(const int m, const int n, const float * d, float * x, const int ldx)
{
	int i, j;

	mc_unused(ldx);
	for (j = 1; j <= n; ++j) {
		for (i = 1; i <= m; ++i) {
			mc_blas_matrix_at(x, ldx, n, i, j) = mc_blas_matrix_at(x, ldx, n, i, j) * mc_blas_vector_at(d, i);
		}
	}
}

#pragma mark - mc_lapack_dlascl2 -

MC_TARGET_FUNC void mc_lapack_dlascl2(const int m, const int n, const double * d, double * x, const int ldx)
{
	int i, j;

	mc_unused(ldx);
	for (j = 1; j <= n; ++j) {
		for (i = 1; i <= m; ++i) {
			mc_blas_matrix_at(x, ldx, n, i, j) = mc_blas_matrix_at(x, ldx, n, i, j) * mc_blas_vector_at(d, i);
		}
	}
}

#pragma mark - mc_lapack_llascl2 -

MC_TARGET_FUNC void mc_lapack_llascl2(const int m, const int n, const long double * d, long double * x, const int ldx)
{
	int i, j;

	mc_unused(ldx);
	for (j = 1; j <= n; ++j) {
		for (i = 1; i <= m; ++i) {
			mc_blas_matrix_at(x, ldx, n, i, j) = mc_blas_matrix_at(x, ldx, n, i, j) * mc_blas_vector_at(d, i);
		}
	}
}

#endif /* !MC_LAPACKE_LASCL2_H */

/* EOF */