//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasq1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/lapack/mc_lapack_lascl.h>
#include <macadam/lapack/mc_lapack_lasq2.h>
#include <macadam/lapack/mc_lapack_lasrt.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_LAPACKE_LASQ1_H
#define MC_LAPACKE_LASQ1_H

#pragma mark - mc_lapack_slasq1 -

MC_TARGET_PROC void mc_lapack_slasq1(const int n, float * d, float * e, float * work, int * info)
{
	const float zero = 0.0f;

	int i, iinfo;
	float eps, scale, safmin, sigmn, sigmx;

	*info = 0;
	if (n < 0) {
		*info = -1;
		mc_blas_xerbla("SLASQ1", -(*info));
		return;
	} else if (n == 0) {
		return;
	} else if (n == 1) {
		mc_blas_vector_at(d, 1) = mc_fabsf(mc_blas_vector_at(d, 1));
		return;
	} else if (n == 2) {
		mc_lapack_slas2(
			  mc_blas_vector_at(d, 1), mc_blas_vector_at(e, 1)
			, mc_blas_vector_at(d, 2)
			, &sigmn
			, &sigmx
		);
		mc_blas_vector_at(d, 1) = sigmx;
		mc_blas_vector_at(d, 2) = sigmn;
		return;
	}

	sigmx = zero;
	for (i = 1; i <= (n - 1); ++i) {
		mc_blas_vector_at(d, i) = mc_fabsf(mc_blas_vector_at(d, i));
		sigmx                   = mc_fmaxf(sigmx, mc_fabsf(mc_blas_vector_at(e, i)));
	}
	mc_blas_vector_at(d, n) = mc_fabsf(mc_blas_vector_at(d, n));

	if (sigmx == zero) {
		mc_lapack_slasrt('D', n, d, &iinfo);
		return;
	}

	for (i = 1; i <= n; ++i) {
		sigmx = mc_fmaxf(sigmx, mc_blas_vector_at(d, i));
	}

	eps    = mc_lapack_slamch('P');
	safmin = mc_lapack_slamch('S');
	scale  = mc_sqrtf(eps / safmin);
	mc_scopy(n, d, 1, work, 2);
	mc_scopy(n - 1, e, 1, work + 1, 2);
	mc_lapack_slascl('G', 0, 0, sigmx, scale, (2 * n) - 1, 1, work, (2 * n) - 1, &iinfo);

	for (i = 1; i <= ((2 * n) - 1); ++i) {
		mc_blas_vector_at(work, i) = mc_raise2f(mc_blas_vector_at(work, i));
	}
	mc_blas_vector_at(work, 2 * n) = zero;
	mc_lapack_slasq2(n, work, info);

	if ((*info) == 0) {
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(d, i) = mc_sqrtf(mc_blas_vector_at(work, i));
		}
		mc_lapack_slascl('G', 0, 0, scale, sigmx, n, 1, d, n, &iinfo);
	} else if ((*info) == 2) {
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(d, i) = mc_sqrtf(mc_blas_vector_at(work, (2 * i) - 1));
			mc_blas_vector_at(e, i) = mc_sqrtf(mc_blas_vector_at(work, i * 2));
		}
		mc_lapack_slascl('G', 0, 0, scale, sigmx, n, 1, d, n, &iinfo);
		mc_lapack_slascl('G', 0, 0, scale, sigmx, n, 1, e, n, &iinfo);
	}
}

#pragma mark - mc_lapack_dlasq1 -

MC_TARGET_PROC void mc_lapack_dlasq1(const int n, double * d, double * e, double * work, int * info)
{
	const double zero = 0.0;

	int i, iinfo;
	double eps, scale, safmin, sigmn, sigmx;

	*info = 0;
	if (n < 0) {
		*info = -1;
		mc_blas_xerbla("DLASQ1", -(*info));
		return;
	} else if (n == 0) {
		return;
	} else if (n == 1) {
		mc_blas_vector_at(d, 1) = mc_fabs(mc_blas_vector_at(d, 1));
		return;
	} else if (n == 2) {
		mc_lapack_dlas2(
			  mc_blas_vector_at(d, 1), mc_blas_vector_at(e, 1)
			, mc_blas_vector_at(d, 2)
			, &sigmn
			, &sigmx
		);
		mc_blas_vector_at(d, 1) = sigmx;
		mc_blas_vector_at(d, 2) = sigmn;
		return;
	}

	sigmx = zero;
	for (i = 1; i <= (n - 1); ++i) {
		mc_blas_vector_at(d, i) = mc_fabs(mc_blas_vector_at(d, i));
		sigmx                   = mc_fmax(sigmx, mc_fabs(mc_blas_vector_at(e, i)));
	}
	mc_blas_vector_at(d, n) = mc_fabs(mc_blas_vector_at(d, n));

	if (sigmx == zero) {
		mc_lapack_dlasrt('D', n, d, &iinfo);
		return;
	}

	for (i = 1; i <= n; ++i) {
		sigmx = mc_fmax(sigmx, mc_blas_vector_at(d, i));
	}

	eps    = mc_lapack_dlamch('P');
	safmin = mc_lapack_dlamch('S');
	scale  = mc_sqrt(eps / safmin);
	mc_dcopy(n, d, 1, work, 2);
	mc_dcopy(n - 1, e, 1, work + 1, 2);
	mc_lapack_dlascl('G', 0, 0, sigmx, scale, (2 * n) - 1, 1, work, (2 * n) - 1, &iinfo);

	for (i = 1; i <= ((2 * n) - 1); ++i) {
		mc_blas_vector_at(work, i) = mc_raise2(mc_blas_vector_at(work, i));
	}
	mc_blas_vector_at(work, 2 * n) = zero;
	mc_lapack_dlasq2(n, work, info);

	if ((*info) == 0) {
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(d, i) = mc_sqrt(mc_blas_vector_at(work, i));
		}
		mc_lapack_dlascl('G', 0, 0, scale, sigmx, n, 1, d, n, &iinfo);
	} else if ((*info) == 2) {
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(d, i) = mc_sqrt(mc_blas_vector_at(work, (2 * i) - 1));
			mc_blas_vector_at(e, i) = mc_sqrt(mc_blas_vector_at(work, i * 2));
		}
		mc_lapack_dlascl('G', 0, 0, scale, sigmx, n, 1, d, n, &iinfo);
		mc_lapack_dlascl('G', 0, 0, scale, sigmx, n, 1, e, n, &iinfo);
	}
}

#pragma mark - mc_lapack_llasq1 -

MC_TARGET_PROC void mc_lapack_llasq1(const int n, long double * d, long double * e, long double * work, int * info)
{
	const long double zero = 0.0L;

	int i, iinfo;
	long double eps, scale, safmin, sigmn, sigmx;

	*info = 0;
	if (n < 0) {
		*info = -1;
		mc_blas_xerbla("LLASQ1", -(*info));
		return;
	} else if (n == 0) {
		return;
	} else if (n == 1) {
		mc_blas_vector_at(d, 1) = mc_fabsl(mc_blas_vector_at(d, 1));
		return;
	} else if (n == 2) {
		mc_lapack_llas2(
			  mc_blas_vector_at(d, 1), mc_blas_vector_at(e, 1)
			, mc_blas_vector_at(d, 2)
			, &sigmn
			, &sigmx
		);
		mc_blas_vector_at(d, 1) = sigmx;
		mc_blas_vector_at(d, 2) = sigmn;
		return;
	}

	sigmx = zero;
	for (i = 1; i <= (n - 1); ++i) {
		mc_blas_vector_at(d, i) = mc_fabsl(mc_blas_vector_at(d, i));
		sigmx                   = mc_fmaxl(sigmx, mc_fabsl(mc_blas_vector_at(e, i)));
	}
	mc_blas_vector_at(d, n) = mc_fabsl(mc_blas_vector_at(d, n));

	if (sigmx == zero) {
		mc_lapack_llasrt('D', n, d, &iinfo);
		return;
	}

	for (i = 1; i <= n; ++i) {
		sigmx = mc_fmaxl(sigmx, mc_blas_vector_at(d, i));
	}

	eps    = mc_lapack_llamch('P');
	safmin = mc_lapack_llamch('S');
	scale  = mc_sqrtl(eps / safmin);
	mc_lcopy(n, d, 1, work, 2);
	mc_lcopy(n - 1, e, 1, work + 1, 2);
	mc_lapack_llascl('G', 0, 0, sigmx, scale, (2 * n) - 1, 1, work, (2 * n) - 1, &iinfo);

	for (i = 1; i <= ((2 * n) - 1); ++i) {
		mc_blas_vector_at(work, i) = mc_raise2l(mc_blas_vector_at(work, i));
	}
	mc_blas_vector_at(work, 2 * n) = zero;
	mc_lapack_llasq2(n, work, info);

	if ((*info) == 0) {
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(d, i) = mc_sqrtl(mc_blas_vector_at(work, i));
		}
		mc_lapack_llascl('G', 0, 0, scale, sigmx, n, 1, d, n, &iinfo);
	} else if ((*info) == 2) {
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(d, i) = mc_sqrtl(mc_blas_vector_at(work, (2 * i) - 1));
			mc_blas_vector_at(e, i) = mc_sqrtl(mc_blas_vector_at(work, i * 2));
		}
		mc_lapack_llascl('G', 0, 0, scale, sigmx, n, 1, d, n, &iinfo);
		mc_lapack_llascl('G', 0, 0, scale, sigmx, n, 1, e, n, &iinfo);
	}
}

#endif /* !MC_LAPACKE_LASQ1_H */

/* EOF */