//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_axpy.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \synopsis
 *    void ?axpy(n, a, x, incx, y, incy)
 *    int           incx, incy, n
 *    real-floating a, x(*), y(*)
 *
 * \purpose
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in the input vectors x and y.
 *    [in]  a    - real-floating. Specifies the scalar alpha.
 *    [in]  x    - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *    [out] y    - real-floating arrays of size at least (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the increment for the elements of `y`. incy must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Linpack.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_abs1.h>
#include <macadam/details/math/mc_cadd.h>
#include <macadam/details/math/mc_cmul.h>

#ifndef MC_BLAS_AXPY_H
#define MC_BLAS_AXPY_H

#pragma mark - mc_blas_saxpy -

MC_TARGET_FUNC void mc_blas_saxpy(const int n, const float a, const float * x, const int incx, float * y, const int incy)
{
	int i, m, ix, iy, mp1;

	if (n <= 0) {
		return;
	}
	if (a == 0.0f) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 4;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (a * mc_blas_vector_at(x, i));
			}
		}
		if (n < 4) {
			return;
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 4) {
			mc_blas_vector_at(y, i    ) = mc_blas_vector_at(y, i    ) + (a * mc_blas_vector_at(x, i    ));
			mc_blas_vector_at(y, i + 1) = mc_blas_vector_at(y, i + 1) + (a * mc_blas_vector_at(x, i + 1));
			mc_blas_vector_at(y, i + 2) = mc_blas_vector_at(y, i + 2) + (a * mc_blas_vector_at(x, i + 2));
			mc_blas_vector_at(y, i + 3) = mc_blas_vector_at(y, i + 3) + (a * mc_blas_vector_at(x, i + 3));
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
			mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (a * mc_blas_vector_at(x, ix));
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_daxpy -

MC_TARGET_FUNC void mc_blas_daxpy(const int n, const double a, const double * x, const int incx, double * y, const int incy)
{
	int i, m, ix, iy, mp1;

	if (n <= 0) {
		return;
	}
	if (a == 0.0) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 4;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (a * mc_blas_vector_at(x, i));
			}
		}
		if (n < 4) {
			return;
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 4) {
			mc_blas_vector_at(y, i    ) = mc_blas_vector_at(y, i    ) + (a * mc_blas_vector_at(x, i    ));
			mc_blas_vector_at(y, i + 1) = mc_blas_vector_at(y, i + 1) + (a * mc_blas_vector_at(x, i + 1));
			mc_blas_vector_at(y, i + 2) = mc_blas_vector_at(y, i + 2) + (a * mc_blas_vector_at(x, i + 2));
			mc_blas_vector_at(y, i + 3) = mc_blas_vector_at(y, i + 3) + (a * mc_blas_vector_at(x, i + 3));
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
			mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (a * mc_blas_vector_at(x, ix));
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_laxpy -

MC_TARGET_FUNC void mc_blas_laxpy(const int n, const long double a, const long double * x, const int incx, long double * y, const int incy)
{
	int i, m, ix, iy, mp1;

	if (n <= 0) {
		return;
	}
	if (a == 0.0L) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 4;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (a * mc_blas_vector_at(x, i));
			}
		}
		if (n < 4) {
			return;
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 4) {
			mc_blas_vector_at(y, i    ) = mc_blas_vector_at(y, i    ) + (a * mc_blas_vector_at(x, i    ));
			mc_blas_vector_at(y, i + 1) = mc_blas_vector_at(y, i + 1) + (a * mc_blas_vector_at(x, i + 1));
			mc_blas_vector_at(y, i + 2) = mc_blas_vector_at(y, i + 2) + (a * mc_blas_vector_at(x, i + 2));
			mc_blas_vector_at(y, i + 3) = mc_blas_vector_at(y, i + 3) + (a * mc_blas_vector_at(x, i + 3));
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
			mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (a * mc_blas_vector_at(x, ix));
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

/* \name
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \synopsis
 *    void ?axpy(n, a, x, incx, y, incy)
 *    int     incx, incy, n
 *    complex a, x(*), y(*)
 *
 * \purpose
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in the input vectors x and y.
 *    [in]  a    - complex. Specifies the scalar alpha.
 *    [in]  x    - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *    [out] y    - complex array of size at least (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the increment for the elements of `y`. incy must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Linpack.
 */

#pragma mark - mc_blas_caxpy -

MC_TARGET_FUNC void mc_blas_caxpy(const int n, const mc_complex_float_t a, const mc_complex_float_t * x, const int incx, mc_complex_float_t * y, const int incy)
{
	int i, m, ix, iy, mp1;

	if (n <= 0) {
		return;
	}
	if (mc_blas_scabs1(a) == 0.0f) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 4;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(y, i) = mc_caddf(mc_blas_vector_at(y, i), mc_cmulf(a, mc_blas_vector_at(x, i)));
			}
		}
		if (n < 4) {
			return;
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 4) {
			mc_blas_vector_at(y, i    ) = mc_caddf(mc_blas_vector_at(y, i    ), mc_cmulf(a, mc_blas_vector_at(x, i    )));
			mc_blas_vector_at(y, i + 1) = mc_caddf(mc_blas_vector_at(y, i + 1), mc_cmulf(a, mc_blas_vector_at(x, i + 1)));
			mc_blas_vector_at(y, i + 2) = mc_caddf(mc_blas_vector_at(y, i + 2), mc_cmulf(a, mc_blas_vector_at(x, i + 2)));
			mc_blas_vector_at(y, i + 3) = mc_caddf(mc_blas_vector_at(y, i + 3), mc_cmulf(a, mc_blas_vector_at(x, i + 3)));
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
			mc_blas_vector_at(y, iy) = mc_caddf(mc_blas_vector_at(y, iy), mc_cmulf(a, mc_blas_vector_at(x, ix)));
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_zaxpy -

MC_TARGET_FUNC void mc_blas_zaxpy(const int n, const mc_complex_double_t a, const mc_complex_double_t * x, const int incx, mc_complex_double_t * y, const int incy)
{
	int i, m, ix, iy, mp1;

	if (n <= 0) {
		return;
	}
	if (mc_blas_dzabs1(a) == 0.0) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 4;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(y, i) = mc_cadd(mc_blas_vector_at(y, i), mc_cmul(a, mc_blas_vector_at(x, i)));
			}
		}
		if (n < 4) {
			return;
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 4) {
			mc_blas_vector_at(y, i    ) = mc_cadd(mc_blas_vector_at(y, i    ), mc_cmul(a, mc_blas_vector_at(x, i    )));
			mc_blas_vector_at(y, i + 1) = mc_cadd(mc_blas_vector_at(y, i + 1), mc_cmul(a, mc_blas_vector_at(x, i + 1)));
			mc_blas_vector_at(y, i + 2) = mc_cadd(mc_blas_vector_at(y, i + 2), mc_cmul(a, mc_blas_vector_at(x, i + 2)));
			mc_blas_vector_at(y, i + 3) = mc_cadd(mc_blas_vector_at(y, i + 3), mc_cmul(a, mc_blas_vector_at(x, i + 3)));
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
			mc_blas_vector_at(y, iy) = mc_cadd(mc_blas_vector_at(y, iy), mc_cmul(a, mc_blas_vector_at(x, ix)));
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_qaxpy -

MC_TARGET_FUNC void mc_blas_qaxpy(const int n, const mc_complex_long_double_t a, const mc_complex_long_double_t * x, const int incx, mc_complex_long_double_t * y, const int incy)
{
	int i, m, ix, iy, mp1;

	if (n <= 0) {
		return;
	}
	if (mc_blas_lqabs1(a) == 0.0L) {
		return;
	}
	if (incx == 1 && incy == 1) {
		m = n % 4;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				mc_blas_vector_at(y, i) = mc_caddl(mc_blas_vector_at(y, i), mc_cmull(a, mc_blas_vector_at(x, i)));
			}
		}
		if (n < 4) {
			return;
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 4) {
			mc_blas_vector_at(y, i    ) = mc_caddl(mc_blas_vector_at(y, i    ), mc_cmull(a, mc_blas_vector_at(x, i    )));
			mc_blas_vector_at(y, i + 1) = mc_caddl(mc_blas_vector_at(y, i + 1), mc_cmull(a, mc_blas_vector_at(x, i + 1)));
			mc_blas_vector_at(y, i + 2) = mc_caddl(mc_blas_vector_at(y, i + 2), mc_cmull(a, mc_blas_vector_at(x, i + 2)));
			mc_blas_vector_at(y, i + 3) = mc_caddl(mc_blas_vector_at(y, i + 3), mc_cmull(a, mc_blas_vector_at(x, i + 3)));
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
			mc_blas_vector_at(y, iy) = mc_caddl(mc_blas_vector_at(y, iy), mc_cmull(a, mc_blas_vector_at(x, ix)));
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#endif /* !MC_BLAS_AXPY_H */

/* EOF */