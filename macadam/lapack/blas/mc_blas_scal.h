//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_scal.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?scal scales a vector by a constant.
 *
 * \synopsis
 *    void ?scal(n, a, x, incx)
 *    int           incx, n
 *    real-floating a
 *    real-floating x(*)
 *
 * \purpose
 *    ?scal scales a vector by a constant.
 *
 * \parameters
 *    [in] n    - int. Specifies the number of elements in the input vector `x`.
 *    [in] a    - real-floating. Specifies the scalar alpha.
 *    [in] x    - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/details/math/mc_cmul.h>

#ifndef MC_BLAS_SCAL_H
#define MC_BLAS_SCAL_H

#pragma mark - mc_blas_sscal -

MC_TARGET_FUNC void mc_blas_sscal(const int n, float a, float * x, const int incx)
{
	int i, m, mp1, nincx;

	if (n <= 0 || incx <= 0) {
		return;
	}
	if (incx == 1) {
		m = n % 5;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(x, i) = a * mc_blas_vector_at(x, i);
			}
			if (n < 5) {
				return;
			}
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 5) {
			mc_blas_vector_at(x, i    ) = a * mc_blas_vector_at(x, i    );
			mc_blas_vector_at(x, i + 1) = a * mc_blas_vector_at(x, i + 1);
			mc_blas_vector_at(x, i + 2) = a * mc_blas_vector_at(x, i + 2);
			mc_blas_vector_at(x, i + 3) = a * mc_blas_vector_at(x, i + 3);
			mc_blas_vector_at(x, i + 4) = a * mc_blas_vector_at(x, i + 4);
		}
	} else {
		nincx = n * incx;
		for (i = 1; incx < 0 ? i >= nincx : i <= nincx; i += incx) {
			mc_blas_vector_at(x, i) = a * mc_blas_vector_at(x, i);
		}
	}
}

#pragma mark - mc_blas_dscal -

MC_TARGET_FUNC void mc_blas_dscal(const int n, double a, double * x, const int incx)
{
	int i, m, mp1, nincx;

	if (n <= 0 || incx <= 0) {
		return;
	}
	if (incx == 1) {
		m = n % 5;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(x, i) = a * mc_blas_vector_at(x, i);
			}
			if (n < 5) {
				return;
			}
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 5) {
			mc_blas_vector_at(x, i    ) = a * mc_blas_vector_at(x, i    );
			mc_blas_vector_at(x, i + 1) = a * mc_blas_vector_at(x, i + 1);
			mc_blas_vector_at(x, i + 2) = a * mc_blas_vector_at(x, i + 2);
			mc_blas_vector_at(x, i + 3) = a * mc_blas_vector_at(x, i + 3);
			mc_blas_vector_at(x, i + 4) = a * mc_blas_vector_at(x, i + 4);
		}
	} else {
		nincx = n * incx;
		for (i = 1; incx < 0 ? i >= nincx : i <= nincx; i += incx) {
			mc_blas_vector_at(x, i) = a * mc_blas_vector_at(x, i);
		}
	}
}

#pragma mark - mc_blas_lscal -

MC_TARGET_FUNC void mc_blas_lscal(const int n, long double a, long double * x, const int incx)
{
	int i, m, mp1, nincx;

	if (n <= 0 || incx <= 0) {
		return;
	}
	if (incx == 1) {
		m = n % 5;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(x, i) = a * mc_blas_vector_at(x, i);
			}
			if (n < 5) {
				return;
			}
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 5) {
			mc_blas_vector_at(x, i    ) = a * mc_blas_vector_at(x, i    );
			mc_blas_vector_at(x, i + 1) = a * mc_blas_vector_at(x, i + 1);
			mc_blas_vector_at(x, i + 2) = a * mc_blas_vector_at(x, i + 2);
			mc_blas_vector_at(x, i + 3) = a * mc_blas_vector_at(x, i + 3);
			mc_blas_vector_at(x, i + 4) = a * mc_blas_vector_at(x, i + 4);
		}
	} else {
		nincx = n * incx;
		for (i = 1; incx < 0 ? i >= nincx : i <= nincx; i += incx) {
			mc_blas_vector_at(x, i) = a * mc_blas_vector_at(x, i);
		}
	}
}

/* \name
 *    ?scal scales a vector by a constant.
 *
 * \synopsis
 *    void ?scal(n, a, x, incx)
 *    int     incx, n
 *    complex a
 *    complex x(*)
 *
 * \purpose
 *    ?scal scales a vector by a constant.
 *
 * \parameters
 *    [in] n    - int. Specifies the number of elements in the input vector `x`.
 *    [in] a    - complex. Specifies the scalar alpha.
 *    [in] x    - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#pragma mark - mc_blas_cscal -

MC_TARGET_FUNC void mc_blas_cscal(const int n, mc_complex_float_t a, mc_complex_float_t * x, const int incx)
{
	int i, nincx;

	if (n <= 0 || incx <= 0) {
		return;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(x, i) = mc_cmulf(a, mc_blas_vector_at(x, i));
		}
	} else {
		nincx = n * incx;
		for (i = 1; incx < 0 ? i >= nincx : i <= nincx; i += incx) {
			mc_blas_vector_at(x, i) = mc_cmulf(a, mc_blas_vector_at(x, i));
		}
	}
}

#pragma mark - mc_blas_zscal -

MC_TARGET_FUNC void mc_blas_zscal(const int n, mc_complex_double_t a, mc_complex_double_t * x, const int incx)
{
	int i, nincx;

	if (n <= 0 || incx <= 0) {
		return;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(x, i) = mc_cmul(a, mc_blas_vector_at(x, i));
		}
	} else {
		nincx = n * incx;
		for (i = 1; incx < 0 ? i >= nincx : i <= nincx; i += incx) {
			mc_blas_vector_at(x, i) = mc_cmul(a, mc_blas_vector_at(x, i));
		}
	}
}

#pragma mark - mc_blas_qscal -

MC_TARGET_FUNC void mc_blas_qscal(const int n, mc_complex_long_double_t a, mc_complex_long_double_t * x, const int incx)
{
	int i, nincx;

	if (n <= 0 || incx <= 0) {
		return;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			mc_blas_vector_at(x, i) = mc_cmull(a, mc_blas_vector_at(x, i));
		}
	} else {
		nincx = n * incx;
		for (i = 1; incx < 0 ? i >= nincx : i <= nincx; i += incx) {
			mc_blas_vector_at(x, i) = mc_cmull(a, mc_blas_vector_at(x, i));
		}
	}
}

#endif /* !MC_BLAS_SCAL_H */

/* EOF */