//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_ger.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?ger performs the rank 1 operation: a=alpha*x*y' + a.
 *
 * \synopsis
 *    void ?ger(m, n, alpha, x, incx, y, incy, a, lda)
 *    real-floating alpha
 *    int            incx, incy, lda, m, n
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?ger performs the rank 1 operation: a=alpha*x*y' + a where alpha is a scalar,
 *    x is an m element vector, y is an n element vector and a is an m by n matrix.
 *
 * \parameters
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in]  x     - real-floating array of dimension (at least) (1+(m-1)*abs(incx)). The incremented
 *    array `x` must contain the m element vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y     - real-floating array of dimension (at least) (1+(n-1)*abs(incy)). The incremented
 *    array `y` must contain the n element vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [out] a     - real-floating array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients. a is overwritten by the updated matrix.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
 *
 * \examples
 *
 * \level 2 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Jeremy Du Croz, Nag Central Office.
 *     \author Sven Hammarling, Nag Central Office.
 *     \author Richard Hanson, Sandia National Labs.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>
#include <macadam/lapack/blas/mc_blas_xerbla.h>
#include <macadam/details/math/mc_cadd.h>
#include <macadam/details/math/mc_ciseq.h>
#include <macadam/details/math/mc_conj.h>
#include <macadam/details/math/mc_cmul.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_GER_H
#define MC_BLAS_GER_H

#pragma mark - mc_blas_sger -

MC_TARGET_FUNC void mc_blas_sger(const int m, const int n, const float alpha, const float * x, const int incx, const float * y, const int incy, float * a, const int lda)
{
	const float zero = 0.0f;

	float temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("SGER  ", info);
		return;
	}

	if (m == 0 || n == 0 || alpha == zero) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (mc_blas_vector_at(y, jy) != zero) {
				temp = alpha * mc_blas_vector_at(y, jy);
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (mc_blas_vector_at(y, jy) != zero) {
				temp = alpha * mc_blas_vector_at(y, jy);
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

#pragma mark - mc_blas_dger -

MC_TARGET_FUNC void mc_blas_dger(const int m, const int n, const double alpha, const double * x, const int incx, const double * y, const int incy, double * a, const int lda)
{
	const double zero = 0.0;

	double temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("DGER  ", info);
		return;
	}

	if (m == 0 || n == 0 || alpha == zero) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (mc_blas_vector_at(y, jy) != zero) {
				temp = alpha * mc_blas_vector_at(y, jy);
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (mc_blas_vector_at(y, jy) != zero) {
				temp = alpha * mc_blas_vector_at(y, jy);
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

#pragma mark - mc_blas_lger -

MC_TARGET_FUNC void mc_blas_lger(const int m, const int n, const long double alpha, const long double * x, const int incx, const long double * y, const int incy, long double * a, const int lda)
{
	const long double zero = 0.0L;

	long double temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("LGER  ", info);
		return;
	}

	if (m == 0 || n == 0 || alpha == zero) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (mc_blas_vector_at(y, jy) != zero) {
				temp = alpha * mc_blas_vector_at(y, jy);
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (mc_blas_vector_at(y, jy) != zero) {
				temp = alpha * mc_blas_vector_at(y, jy);
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

/* \name
 *    ?gerc performs the rank 1 operation: a=alpha*x*y_ + a.
 *
 * \synopsis
 *    void ?gerc(m, n, alpha, x, incx, y, incy, a, lda)
 *    complex alpha
 *    int     incx, incy, lda, m, n
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?ger performs the rank 1 operation: a=alpha*x*y_ + a where alpha is a scalar,
 *    x is an m element vector, y is an n element vector and a is an m by n matrix.
 *
 * \parameters
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - complex. Specifies the scalar alpha.
 *
 *    [in]  x     - complex array of dimension (at least) (1+(m-1)*abs(incx)). The incremented
 *    array `x` must contain the m element vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y     - complex array of dimension (at least) (1+(n-1)*abs(incy)). The incremented
 *    array `y` must contain the n element vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [out] a     - complex array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients. a is overwritten by the updated matrix.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
 *
 * \examples
 *
 * \level 2 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Jeremy Du Croz, Nag Central Office.
 *     \author Sven Hammarling, Nag Central Office.
 *     \author Richard Hanson, Sandia National Labs.
 */

#pragma mark - mc_blas_cgerc -

MC_TARGET_FUNC void mc_blas_cgerc(const int m, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y, const int incy, mc_complex_float_t * a, const int lda)
{
	const mc_complex_float_t zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("CGERC ", info);
		return;
	}

	if (m == 0 || n == 0 || mc_ciseqf(alpha, zero)) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseqf(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmulf(alpha, mc_conjf(mc_blas_vector_at(y, jy)));
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, i), temp));
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseqf(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmulf(alpha, mc_conjf(mc_blas_vector_at(y, jy)));
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, ix), temp));
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

#pragma mark - mc_blas_zgerc -

MC_TARGET_FUNC void mc_blas_zgerc(const int m, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y, const int incy, mc_complex_double_t * a, const int lda)
{
	const mc_complex_double_t zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("ZGERC ", info);
		return;
	}

	if (m == 0 || n == 0 || mc_ciseq(alpha, zero)) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseq(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmul(alpha, mc_conj(mc_blas_vector_at(y, jy)));
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, i), temp));
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseq(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmul(alpha, mc_conj(mc_blas_vector_at(y, jy)));
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, ix), temp));
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

#pragma mark - mc_blas_qgerc -

MC_TARGET_FUNC void mc_blas_qgerc(const int m, const int n, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * x, const int incx, const mc_complex_long_double_t * y, const int incy, mc_complex_long_double_t * a, const int lda)
{
	const mc_complex_long_double_t zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("QGERC ", info);
		return;
	}

	if (m == 0 || n == 0 || mc_ciseql(alpha, zero)) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseql(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmull(alpha, mc_conjl(mc_blas_vector_at(y, jy)));
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, i), temp));
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseql(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmull(alpha, mc_conjl(mc_blas_vector_at(y, jy)));
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, ix), temp));
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

/* \name
 *    ?geru performs the rank 1 operation: a=alpha*x*y' + a.
 *
 * \synopsis
 *    void ?geru(m, n, alpha, x, incx, y, incy, a, lda)
 *    complex alpha
 *    int     incx, incy, lda, m, n
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?ger performs the rank 1 operation: a=alpha*x*y' + a where alpha is a scalar,
 *    x is an m element vector, y is an n element vector and a is an m by n matrix.
 *
 * \parameters
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - complex. Specifies the scalar alpha.
 *
 *    [in]  x     - complex array of dimension (at least) (1+(m-1)*abs(incx)). The incremented
 *    array `x` must contain the m element vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y     - complex array of dimension (at least) (1+(n-1)*abs(incy)). The incremented
 *    array `y` must contain the n element vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [out] a     - complex array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients. a is overwritten by the updated matrix.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
 *
 * \examples
 *
 * \level 2 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Jeremy Du Croz, Nag Central Office.
 *     \author Sven Hammarling, Nag Central Office.
 *     \author Richard Hanson, Sandia National Labs.
 */

#pragma mark - mc_blas_cgeru -

MC_TARGET_FUNC void mc_blas_cgeru(const int m, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y, const int incy, mc_complex_float_t * a, const int lda)
{
	const mc_complex_float_t zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("CGERU ", info);
		return;
	}

	if (m == 0 || n == 0 || mc_ciseqf(alpha, zero)) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseqf(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmulf(alpha, mc_blas_vector_at(y, jy));
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, i), temp));
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseqf(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmulf(alpha, mc_blas_vector_at(y, jy));
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, ix), temp));
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

#pragma mark - mc_blas_zgeru -

MC_TARGET_FUNC void mc_blas_zgeru(const int m, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y, const int incy, mc_complex_double_t * a, const int lda)
{
	const mc_complex_double_t zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("ZGERU ", info);
		return;
	}

	if (m == 0 || n == 0 || mc_ciseq(alpha, zero)) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseq(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmul(alpha, mc_blas_vector_at(y, jy));
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, i), temp));
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseq(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmul(alpha, mc_blas_vector_at(y, jy));
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, ix), temp));
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

#pragma mark - mc_blas_qgeru -

MC_TARGET_FUNC void mc_blas_qgeru(const int m, const int n, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * x, const int incx, const mc_complex_long_double_t * y, const int incy, mc_complex_long_double_t * a, const int lda)
{
	const mc_complex_long_double_t zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp;
	int i, info, j, ix, jy, kx;

	info = 0;
	if (m < 0) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("QGERU ", info);
		return;
	}

	if (m == 0 || n == 0 || mc_ciseql(alpha, zero)) {
		return;
	}

	if (incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (n - 1) * incy;
	}
	if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseql(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmull(alpha, mc_blas_vector_at(y, jy));
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, i), temp));
				}
			}
			jy = jy + incy;
		}
	} else {
		if (incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (m - 1) * incx;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			if (!mc_ciseql(mc_blas_vector_at(y, jy), zero)) {
				temp = mc_cmull(alpha, mc_blas_vector_at(y, jy));
				ix   = kx;
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, ix), temp));
					ix                                 = ix + incx;
				}
			}
			jy = jy + incy;
		}
	}
}

#endif /* !MC_BLAS_GER_H */

/* EOF */