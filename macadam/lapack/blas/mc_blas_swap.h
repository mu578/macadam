//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_swap.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?swap interchanges two vectors.
 *
 * \synopsis
 *    void ?swap(n, x, incx, y, incy)
 *    int           incx, incy, n
 *    real-floating x(*), y(*)
 *
 * \purpose
 *    ?swap interchanges two vectors.
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in input vector(s).
 *
 *    [out] x    - real-floating array of dimension (at least) (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the storage spacing between elements of `x`.
 *
 *    [out]  y   - real-floating array of dimension (at least) (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the storage spacing between elements of `y`.
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

#ifndef MC_BLAS_SWAP_H
#define MC_BLAS_SWAP_H

#pragma mark - mc_blas_sswap -

MC_TARGET_FUNC void mc_blas_sswap(const int n, float * x, const int incx, float * y, const int incy)
{
	float temp;
	int i, ix, iy, m, mp1;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 3;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				temp                    = mc_blas_vector_at(x, i);
				mc_blas_vector_at(x, i) = mc_blas_vector_at(y, i);
				mc_blas_vector_at(y, i) = temp;
			}
			if (n < 3) {
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
		for (i = mp1; i <= n; i += 3) {
			temp                        = mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i    ) = mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i    ) = temp;
			temp                        = mc_blas_vector_at(x, i + 1);
			mc_blas_vector_at(x, i + 1) = mc_blas_vector_at(y, i + 1);
			mc_blas_vector_at(y, i + 1) = temp;
			temp                        = mc_blas_vector_at(x, i + 2);
			mc_blas_vector_at(x, i + 2) = mc_blas_vector_at(y, i + 2);
			mc_blas_vector_at(y, i + 2) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_dswap -

MC_TARGET_FUNC void mc_blas_dswap(const int n, double * x, const int incx, double * y, const int incy)
{
	double temp;
	int i, ix, iy, m, mp1;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 3;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				temp                    = mc_blas_vector_at(x, i);
				mc_blas_vector_at(x, i) = mc_blas_vector_at(y, i);
				mc_blas_vector_at(y, i) = temp;
			}
			if (n < 3) {
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
		for (i = mp1; i <= n; i += 3) {
			temp                        = mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i    ) = mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i    ) = temp;
			temp                        = mc_blas_vector_at(x, i + 1);
			mc_blas_vector_at(x, i + 1) = mc_blas_vector_at(y, i + 1);
			mc_blas_vector_at(y, i + 1) = temp;
			temp                        = mc_blas_vector_at(x, i + 2);
			mc_blas_vector_at(x, i + 2) = mc_blas_vector_at(y, i + 2);
			mc_blas_vector_at(y, i + 2) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_lswap -

MC_TARGET_FUNC void mc_blas_lswap(const int n, long double * x, const int incx, long double * y, const int incy)
{
	long double temp;
	int i, ix, iy, m, mp1;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 3;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				temp                    = mc_blas_vector_at(x, i);
				mc_blas_vector_at(x, i) = mc_blas_vector_at(y, i);
				mc_blas_vector_at(y, i) = temp;
			}
			if (n < 3) {
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
		for (i = mp1; i <= n; i += 3) {
			temp                        = mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i    ) = mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i    ) = temp;
			temp                        = mc_blas_vector_at(x, i + 1);
			mc_blas_vector_at(x, i + 1) = mc_blas_vector_at(y, i + 1);
			mc_blas_vector_at(y, i + 1) = temp;
			temp                        = mc_blas_vector_at(x, i + 2);
			mc_blas_vector_at(x, i + 2) = mc_blas_vector_at(y, i + 2);
			mc_blas_vector_at(y, i + 2) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

/* \name
 *    ?swap interchanges two vectors.
 *
 * \synopsis
 *    void ?swap(n, x, incx, y, incy)
 *    int     incx, incy, n
 *    complex x(*), y(*)
 *
 * \purpose
 *    ?swap interchanges two vectors.
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in input vector(s).
 *
 *    [out] x    - complex array of dimension (at least) (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the storage spacing between elements of `x`.
 *
 *    [out]  y   - complex array of dimension (at least) (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the storage spacing between elements of `y`.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#pragma mark - mc_blas_cswap -

MC_TARGET_FUNC void mc_blas_cswap(const int n, mc_complex_float_t * x, const int incx, mc_complex_float_t * y, const int incy)
{
	mc_complex_float_t temp;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i) = mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_zswap -

MC_TARGET_FUNC void mc_blas_zswap(const int n, mc_complex_double_t * x, const int incx, mc_complex_double_t * y, const int incy)
{
	mc_complex_double_t temp;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i) = mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_qswap -

MC_TARGET_FUNC void mc_blas_qswap(const int n, mc_complex_long_double_t * x, const int incx, mc_complex_long_double_t * y, const int incy)
{
	mc_complex_long_double_t temp;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i) = mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#endif /* !MC_BLAS_SWAP_H */

/* EOF */