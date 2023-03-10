//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_syr.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?syr performs the symmetric rank 1 operation:
 *    a=alpha*x*x' + a.
 *
 * \synopsis
 *    void ?syr(uplo, n, alpha, x, incx, a, lda)
 *    real-floating alpha
 *    int           incx, lda, n
 *    char          uplo
 *    real-floating a(lda, *), x(*)
 *
 * \purpose
 *    ?syr performs the symmetric rank 1 operation: a=alpha*x*x' + a where alpha is a scalar,
 *    `x` is an n element vector and `a` is an n by n symmetric matrix.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `a`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `a` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `a` is to be referenced.
 *
 *    [in] n     - int. Specifies the order of the matrix `a`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in]  x    - real-floating array of dimension (at least) (1+(n-1)*abs(incx)). The incremented array `x`
 *    must contain the vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [out] a    - real-floating array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper
 *    triangular part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced.
 *    The upper triangular part of the array `a` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower
 *    triangular part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *    The lower triangular part of the array `a` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least max(1, n).
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
#include <macadam/details/math/mc_cmul.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_SYR_H
#define MC_BLAS_SYR_H

#pragma mark - mc_blas_ssyr -

MC_TARGET_FUNC void mc_blas_ssyr(const char uplo, const int n, const float alpha, const float * x, const int incx, float * a, const int lda)
{
	const float zero = 0.0f;

	float temp;
	int i, info, ix, j, jx, kx;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, n)) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("SSYR  ", info);
		return;
	}

	if (n == 0 || alpha == zero) {
		return;
	}

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero) {
					temp = alpha * mc_blas_vector_at(x, j);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, jx) != zero) {
					temp = alpha * mc_blas_vector_at(x, jx);
					ix   = kx;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	} else {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero) {
					temp = alpha * mc_blas_vector_at(x, j);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, jx) != zero) {
					temp = alpha * mc_blas_vector_at(x, jx);
					ix   = jx;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	}
}

#pragma mark - mc_blas_dsyr -

MC_TARGET_FUNC void mc_blas_dsyr(const char uplo, const int n, const double alpha, const double * x, const int incx, double * a, const int lda)
{
	const double zero = 0.0;

	double temp;
	int i, info, ix, j, jx, kx;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, n)) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("DSYR  ", info);
		return;
	}

	if (n == 0 || alpha == zero) {
		return;
	}

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero) {
					temp = alpha * mc_blas_vector_at(x, j);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, jx) != zero) {
					temp = alpha * mc_blas_vector_at(x, jx);
					ix   = kx;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	} else {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero) {
					temp = alpha * mc_blas_vector_at(x, j);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, jx) != zero) {
					temp = alpha * mc_blas_vector_at(x, jx);
					ix   = jx;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	}
}

#pragma mark - mc_blas_lsyr -

MC_TARGET_FUNC void mc_blas_lsyr(const char uplo, const int n, const long double alpha, const long double * x, const int incx, long double * a, const int lda)
{
	const long double zero = 0.0L;

	long double temp;
	int i, info, ix, j, jx, kx;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, n)) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("LSYR  ", info);
		return;
	}

	if (n == 0 || alpha == zero) {
		return;
	}

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero) {
					temp = alpha * mc_blas_vector_at(x, j);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, jx) != zero) {
					temp = alpha * mc_blas_vector_at(x, jx);
					ix   = kx;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	} else {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero) {
					temp = alpha * mc_blas_vector_at(x, j);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, i) * temp);
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, jx) != zero) {
					temp = alpha * mc_blas_vector_at(x, jx);
					ix   = jx;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + (mc_blas_vector_at(x, ix) * temp);
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	}
}

/* \name
 *    ?syr performs the symmetric rank 1 operation:
 *    a=alpha*x*x_ + a.
 *
 * \synopsis
 *    void ?syr(uplo, n, alpha, x, incx, a, lda)
 *    complex alpha
 *    int     incx, lda, n
 *    char    uplo
 *    complex a(lda, *), x(*)
 *
 * \purpose
 *    ?syr performs the symmetric rank 1 operation: a=alpha*x*x_ + a where alpha is a scalar,
 *    `x` is an n element vector and `a` is an n by n symmetric matrix.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `a`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `a` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `a` is to be referenced.
 *
 *    [in] n     - int. Specifies the order of the matrix `a`, n must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in]  x    - complex array of dimension (at least) (1+(n-1)*abs(incx)). The incremented array `x`
 *    must contain the vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [out] a    - complex array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper
 *    triangular part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced.
 *    The upper triangular part of the array `a` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower
 *    triangular part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *    The lower triangular part of the array `a` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least max(1, n).
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

#pragma mark - mc_blas_csyr -

MC_TARGET_FUNC void mc_blas_csyr(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, mc_complex_float_t * a, const int lda)
{
	const mc_complex_float_t zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp;
	int i, info, ix, j, jx, kx;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, n)) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("CSYR  ", info);
		return;
	}

	if (n == 0 || mc_ciseqf(alpha, zero)) {
		return;
	}

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, j);
				if (!mc_ciseqf(temp, zero)) {
					temp = mc_cmulf(alpha, temp);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, i), temp));
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, jx);
				if (!mc_ciseqf(temp, zero)) {
					temp = mc_cmulf(alpha, temp);
					ix   = kx;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, ix), temp));
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	} else {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, j);
				if (!mc_ciseqf(temp, zero)) {
					temp = mc_cmulf(alpha, temp);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, i), temp));
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, jx);
				if (!mc_ciseqf(temp, zero)) {
					temp = mc_cmulf(alpha, temp);
					ix   = jx;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_cmulf(mc_blas_vector_at(x, ix), temp));
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	}
}

#pragma mark - mc_blas_zsyr -

MC_TARGET_FUNC void mc_blas_zsyr(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, mc_complex_double_t * a, const int lda)
{
	const mc_complex_double_t zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp;
	int i, info, ix, j, jx, kx;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, n)) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("ZSYR  ", info);
		return;
	}

	if (n == 0 || mc_ciseq(alpha, zero)) {
		return;
	}

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, j);
				if (!mc_ciseq(temp, zero)) {
					temp = mc_cmul(alpha, temp);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, i), temp));
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, jx);
				if (!mc_ciseq(temp, zero)) {
					temp = mc_cmul(alpha, temp);
					ix   = kx;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, ix), temp));
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	} else {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, j);
				if (!mc_ciseq(temp, zero)) {
					temp = mc_cmul(alpha, temp);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, i), temp));
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, jx);
				if (!mc_ciseq(temp, zero)) {
					temp = mc_cmul(alpha, temp);
					ix   = jx;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cmul(mc_blas_vector_at(x, ix), temp));
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	}
}

#pragma mark - mc_blas_qsyr -

MC_TARGET_FUNC void mc_blas_qsyr(const char uplo, const int n, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * x, const int incx, mc_complex_long_double_t * a, const int lda)
{
	const mc_complex_long_double_t zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp;
	int i, info, ix, j, jx, kx;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, n)) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("QSYR  ", info);
		return;
	}

	if (n == 0 || mc_ciseql(alpha, zero)) {
		return;
	}

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, j);
				if (!mc_ciseql(temp, zero)) {
					temp = mc_cmull(alpha, temp);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, i), temp));
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, jx);
				if (!mc_ciseql(temp, zero)) {
					temp = mc_cmull(alpha, temp);
					ix   = kx;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, ix), temp));
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	} else {
		if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, j);
				if (!mc_ciseql(temp, zero)) {
					temp = mc_cmull(alpha, temp);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, i), temp));
					}
				}
			}
		} else {
			jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp = mc_blas_vector_at(x, jx);
				if (!mc_ciseql(temp, zero)) {
					temp = mc_cmull(alpha, temp);
					ix   = jx;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_cmull(mc_blas_vector_at(x, ix), temp));
						ix                                 = ix + incx;
					}
				}
				jx = jx + incx;
			}
		}
	}
}

#endif /* !MC_BLAS_SYR_H */

/* EOF */