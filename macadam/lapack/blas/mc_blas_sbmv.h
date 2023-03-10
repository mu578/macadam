//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_sbmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?sbmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?sbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
 *    real-floating alpha, beta
 *    int           incx, incy, k, lda, n
 *    char          uplo
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?sbmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `a` is an n by n symmetric band matrix, with k
 *    super-diagonals. It computes the matrix-vector product for a real symmetric band matrix,
 *    The band matrix A is stored in either upper or lower-band-packed storage mode, it uses
 *    the scalars alpha and beta, vectors x and y, and band matrix `a`.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the band matrix `a`
 *    is being supplied as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] n     - int. Specifies the order of the symmetric matrix `a`, n must be at least zero.
 *    [in] k     - int. Specifies the number of super-diagonals of the matrix symmetric `a`, k
 *    must satisfy 0 < k, i.e must be at least one.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row (k + 1) of the array, the first super-diagonal starting at position 1 in row k, and so on. The
 *    top left k by k triangle of the array `a` is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the lower triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row 0 of the array, the first sub-diagonal starting at position 0 in row 1, and so on. The bottom
 *    right k by k triangle of the array `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. `a` must be at least (k + 1).
 *
 *    [int] x    - real-floating array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must
 *    contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - real-floating. Specifies the scalar beta.
 *
 *    [out] y    - real-floating array of size at least (1+(n-1)*abs(incy)). The incremented array `y` must
 *    contain the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 * \examples
 *              | 1 | 1 | 1 | 1 | 0 | 0 | 0 |
 *              | 1 | 2 | 2 | 2 | 2 | 0 | 0 |
 *              | 1 | 2 | 3 | 3 | 3 | 3 | 0 |
 *     a[7x7] = | 1 | 2 | 3 | 4 | 4 | 4 | 4 |
 *              | 0 | 2 | 3 | 4 | 5 | 5 | 5 |
 *              | 0 | 0 | 3 | 4 | 5 | 6 | 6 |
 *              | 0 | 0 | 0 | 4 | 5 | 6 | 7 |
 *
 *     const real-floating a_band[] = {
 *          0, 0, 0, 1, 2, 3, 4
 *        , 0, 0, 1, 2, 3, 4, 5
 *        , 0, 1, 2, 3, 4, 5, 6
 *        , 1, 2, 3, 4, 5, 6, 7
 *        , 0, 0, 0, 0, 0, 0, 0
 *     };
 *     const real-floating x[] = { 1, 2, 3, 4, 5, 6, 7 };
 *           real-floating y[] = { 1, 0 , 2, 0 , 3, 0 , 4, 0 , 5, 0 , 6, 0 , 7 };
 *     mc_blas_?sbmv('U', 7, 3, 2, a_band, 5, x, 1, 10, y, 2);
 *     on output -> y = { 30, 0 , 78, 0 , 148, 0 , 244, 0 , 288, 0 , 316, 0 , 322 }
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
#include <macadam/details/math/mc_minmag.h>

#ifndef MC_BLAS_SBMV_H
#define MC_BLAS_SBMV_H

#pragma mark - mc_blas_ssbmv -

MC_TARGET_FUNC void mc_blas_ssbmv(const char uplo, const int n, const int k, const float alpha, const float * a, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)
{
	const float one = 1.0f, zero = 0.0f;

	float temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (k < 0) {
		info = 3;
	} else if (lda < k + 1) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	} else if (incy == 0) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("SSBMV ", info);
		return;
	}

	if (n == 0 || (alpha == zero && beta == one)) {
		return;
	}

	if (incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (n - 1) * incx;
	}
	if (incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (n - 1) * incy;
	}

	if (beta != one) {
		if (incy == 1) {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = zero;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = beta * mc_blas_vector_at(y, i);
				}
			}
		} else {
			iy = ky;
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = zero;
					iy                       = iy + incy;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = beta * mc_blas_vector_at(y, iy);
					iy                       = iy + incy;
				}
			}
		}
	}

	if (alpha == zero) {
		return;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		kplus1 = k + 1;
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = alpha * mc_blas_vector_at(x, j);
				temp2 = zero;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                   = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + temp1 * mc_blas_matrix_at(a, lda, n, kplus1, j) + alpha * temp2;
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = alpha * mc_blas_vector_at(x, jx);
				temp2 = zero;
				ix    = kx;
				iy    = ky;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                    = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + temp1 * mc_blas_matrix_at(a, lda, n, kplus1, j) + alpha * temp2;
				jx                       = jx + incx;
				jy                       = jy + incy;
				if (j > k) {
					kx = kx + incx;
					ky = ky + incy;
				}
			}
		}
	} else {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                   = alpha * mc_blas_vector_at(x, j);
				temp2                   = zero;
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (temp1 * mc_blas_matrix_at(a, lda, n, 1, j));
				l                       = 1 - j;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                   = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (alpha * temp2);
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                    = alpha * mc_blas_vector_at(x, jx);
				temp2                    = zero;
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (temp1 * mc_blas_matrix_at(a, lda, n, 1, j));
				l                        = 1 - j;
				ix                       = jx;
				iy                       = jy;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                    = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (alpha * temp2);
				jx                       = jx + incx;
				jy                       = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_dsbmv -

MC_TARGET_FUNC void mc_blas_dsbmv(const char uplo, const int n, const int k, const double alpha, const double * a, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)
{
	const double one = 1.0, zero = 0.0;

	double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (k < 0) {
		info = 3;
	} else if (lda < k + 1) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	} else if (incy == 0) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("DSBMV ", info);
		return;
	}

	if (n == 0 || (alpha == zero && beta == one)) {
		return;
	}

	if (incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (n - 1) * incx;
	}
	if (incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (n - 1) * incy;
	}

	if (beta != one) {
		if (incy == 1) {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = zero;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = beta * mc_blas_vector_at(y, i);
				}
			}
		} else {
			iy = ky;
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = zero;
					iy                       = iy + incy;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = beta * mc_blas_vector_at(y, iy);
					iy                       = iy + incy;
				}
			}
		}
	}

	if (alpha == zero) {
		return;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		kplus1 = k + 1;
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = alpha * mc_blas_vector_at(x, j);
				temp2 = zero;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                   = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + temp1 * mc_blas_matrix_at(a, lda, n, kplus1, j) + alpha * temp2;
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = alpha * mc_blas_vector_at(x, jx);
				temp2 = zero;
				ix    = kx;
				iy    = ky;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                    = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + temp1 * mc_blas_matrix_at(a, lda, n, kplus1, j) + alpha * temp2;
				jx                       = jx + incx;
				jy                       = jy + incy;
				if (j > k) {
					kx = kx + incx;
					ky = ky + incy;
				}
			}
		}
	} else {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                   = alpha * mc_blas_vector_at(x, j);
				temp2                   = zero;
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (temp1 * mc_blas_matrix_at(a, lda, n, 1, j));
				l                       = 1 - j;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                   = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (alpha * temp2);
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                    = alpha * mc_blas_vector_at(x, jx);
				temp2                    = zero;
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (temp1 * mc_blas_matrix_at(a, lda, n, 1, j));
				l                        = 1 - j;
				ix                       = jx;
				iy                       = jy;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                    = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (alpha * temp2);
				jx                       = jx + incx;
				jy                       = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_lsbmv -

MC_TARGET_FUNC void mc_blas_lsbmv(const char uplo, const int n, const int k, const long double alpha, const long double * a, const int lda, const long double * x, const int incx, const long double beta, long double * y, const int incy)
{
	const long double one = 1.0, zero = 0.0;

	long double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (k < 0) {
		info = 3;
	} else if (lda < k + 1) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	} else if (incy == 0) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("LSBMV ", info);
		return;
	}

	if (n == 0 || (alpha == zero && beta == one)) {
		return;
	}

	if (incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (n - 1) * incx;
	}
	if (incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (n - 1) * incy;
	}

	if (beta != one) {
		if (incy == 1) {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = zero;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = beta * mc_blas_vector_at(y, i);
				}
			}
		} else {
			iy = ky;
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = zero;
					iy                       = iy + incy;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = beta * mc_blas_vector_at(y, iy);
					iy                       = iy + incy;
				}
			}
		}
	}

	if (alpha == zero) {
		return;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		kplus1 = k + 1;
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = alpha * mc_blas_vector_at(x, j);
				temp2 = zero;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                   = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + temp1 * mc_blas_matrix_at(a, lda, n, kplus1, j) + alpha * temp2;
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = alpha * mc_blas_vector_at(x, jx);
				temp2 = zero;
				ix    = kx;
				iy    = ky;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                    = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + temp1 * mc_blas_matrix_at(a, lda, n, kplus1, j) + alpha * temp2;
				jx                       = jx + incx;
				jy                       = jy + incy;
				if (j > k) {
					kx = kx + incx;
					ky = ky + incy;
				}
			}
		}
	} else {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                   = alpha * mc_blas_vector_at(x, j);
				temp2                   = zero;
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (temp1 * mc_blas_matrix_at(a, lda, n, 1, j));
				l                       = 1 - j;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                   = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (alpha * temp2);
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                    = alpha * mc_blas_vector_at(x, jx);
				temp2                    = zero;
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (temp1 * mc_blas_matrix_at(a, lda, n, 1, j));
				l                        = 1 - j;
				ix                       = jx;
				iy                       = jy;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_matrix_at(a, lda, n, l + i, j));
					temp2                    = temp2 + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (alpha * temp2);
				jx                       = jx + incx;
				jy                       = jy + incy;
			}
		}
	}
}

/* \name
 *    ?sbmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?sbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
 *    complex alpha, beta
 *    int     incx, incy, k, lda, n
 *    char    uplo
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?sbmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `a` is an n by n symmetric band matrix, with k
 *    super-diagonals. It computes the matrix-vector product for a real symmetric band matrix,
 *    The band matrix A is stored in either upper or lower-band-packed storage mode, it uses
 *    the scalars alpha and beta, vectors x and y, and band matrix `a`.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the band matrix `a`
 *    is being supplied as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] n     - int. Specifies the order of the symmetric matrix `a`, n must be at least zero.
 *    [in] k     - int. Specifies the number of super-diagonals of the matrix symmetric `a`, k
 *    must satisfy 0 < k, i.e must be at least one.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row (k + 1) of the array, the first super-diagonal starting at position 1 in row k, and so on. The
 *    top left k by k triangle of the array `a` is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the lower triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row 0 of the array, the first sub-diagonal starting at position 0 in row 1, and so on. The bottom
 *    right k by k triangle of the array `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. `a` must be at least (k + 1).
 *
 *    [int] x    - complex array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must
 *    contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - complex. Specifies the scalar beta.
 *
 *    [out] y    - complex array of size at least (1+(n-1)*abs(incy)). The incremented array `y` must
 *    contain the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
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

#pragma mark - mc_blas_csbmv -

MC_TARGET_FUNC void mc_blas_csbmv(const char uplo, const int n, const int k, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * x, const int incx, const mc_complex_float_t beta, mc_complex_float_t * y, const int incy)
{
	const mc_complex_float_t one = mc_cmplxf(1.0f, 0.0f), zero =  mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (k < 0) {
		info = 3;
	} else if (lda < k + 1) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	} else if (incy == 0) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("CSBMV ", info);
		return;
	}

	if (n == 0 || (mc_ciseqf(alpha, zero) && mc_ciseqf(beta, one))) {
		return;
	}

	if (incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (n - 1) * incx;
	}
	if (incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (n - 1) * incy;
	}

	if (!mc_ciseqf(beta, one)) {
		if (incy == 1) {
			if (mc_ciseqf(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = zero;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_cmulf(beta, mc_blas_vector_at(y, i));
				}
			}
		} else {
			iy = ky;
			if (mc_ciseqf(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = zero;
					iy                       = iy + incy;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = mc_cmulf(beta, mc_blas_vector_at(y, iy));
					iy                       = iy + incy;
				}
			}
		}
	}

	if (mc_ciseqf(alpha, zero)) {
		return;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		kplus1 = k + 1;
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_cmulf(alpha, mc_blas_vector_at(x, j));
				temp2 = zero;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_caddf(mc_blas_vector_at(y, i), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                   = mc_caddf(temp2, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
				}
				mc_blas_vector_at(y, j) = mc_caddf(mc_blas_vector_at(y, j), mc_caddf(mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, kplus1, j)), mc_cmulf(alpha, temp2)));
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_cmulf(alpha, mc_blas_vector_at(x, jx));
				temp2 = zero;
				ix    = kx;
				iy    = ky;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, iy) = mc_caddf(mc_blas_vector_at(y, iy), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                    = mc_caddf(temp2, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_caddf(mc_blas_vector_at(y, jy), mc_caddf(mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, kplus1, j)), mc_cmulf(alpha, temp2)));
				jx                       = jx + incx;
				jy                       = jy + incy;
				if (j > k) {
					kx = kx + incx;
					ky = ky + incy;
				}
			}
		}
	} else {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                   = mc_cmulf(alpha, mc_blas_vector_at(x, j));
				temp2                   = zero;
				mc_blas_vector_at(y, j) = mc_caddf(mc_blas_vector_at(y, j), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, 1, j)));
				l                       = 1 - j;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					mc_blas_vector_at(y, i) = mc_caddf(mc_blas_vector_at(y, i), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                   = mc_caddf(temp2, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
				}
				mc_blas_vector_at(y, j) = mc_caddf(mc_blas_vector_at(y, j), mc_cmulf(alpha, temp2));
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                    = mc_cmulf(alpha, mc_blas_vector_at(x, jx));
				temp2                    = zero;
				mc_blas_vector_at(y, jy) = mc_caddf(mc_blas_vector_at(y, jy), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, 1, j)));
				l                        = 1 - j;
				ix                       = jx;
				iy                       = jy;
				for (i = (j + 1); i <= mc_minmag(n , j + k); ++i) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_caddf(mc_blas_vector_at(y, iy), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                    = mc_caddf(temp2, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
				}
				mc_blas_vector_at(y, jy) = mc_caddf(mc_blas_vector_at(y, jy), mc_cmulf(alpha, temp2));
				jx                       = jx + incx;
				jy                       = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_zsbmv -

MC_TARGET_FUNC void mc_blas_zsbmv(const char uplo, const int n, const int k, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * x, const int incx, const mc_complex_double_t beta, mc_complex_double_t * y, const int incy)
{
	const mc_complex_double_t one = mc_cmplx(1.0, 0.0), zero =  mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (k < 0) {
		info = 3;
	} else if (lda < k + 1) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	} else if (incy == 0) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("ZSBMV ", info);
		return;
	}

	if (n == 0 || (mc_ciseq(alpha, zero) && mc_ciseq(beta, one))) {
		return;
	}

	if (incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (n - 1) * incx;
	}
	if (incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (n - 1) * incy;
	}

	if (!mc_ciseq(beta, one)) {
		if (incy == 1) {
			if (mc_ciseq(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = zero;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_cmul(beta, mc_blas_vector_at(y, i));
				}
			}
		} else {
			iy = ky;
			if (mc_ciseq(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = zero;
					iy                       = iy + incy;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = mc_cmul(beta, mc_blas_vector_at(y, iy));
					iy                       = iy + incy;
				}
			}
		}
	}

	if (mc_ciseq(alpha, zero)) {
		return;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		kplus1 = k + 1;
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_cmul(alpha, mc_blas_vector_at(x, j));
				temp2 = zero;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_cadd(mc_blas_vector_at(y, i), mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                   = mc_cadd(temp2, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
				}
				mc_blas_vector_at(y, j) = mc_cadd(mc_blas_vector_at(y, j), mc_cadd(mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, kplus1, j)), mc_cmul(alpha, temp2)));
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_cmul(alpha, mc_blas_vector_at(x, jx));
				temp2 = zero;
				ix    = kx;
				iy    = ky;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, iy) = mc_cadd(mc_blas_vector_at(y, iy), mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                    = mc_cadd(temp2, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_cadd(mc_blas_vector_at(y, jy), mc_cadd(mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, kplus1, j)), mc_cmul(alpha, temp2)));
				jx                       = jx + incx;
				jy                       = jy + incy;
				if (j > k) {
					kx = kx + incx;
					ky = ky + incy;
				}
			}
		}
	} else {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                   = mc_cmul(alpha, mc_blas_vector_at(x, j));
				temp2                   = zero;
				mc_blas_vector_at(y, j) = mc_cadd(mc_blas_vector_at(y, j), mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, 1, j)));
				l                       = 1 - j;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					mc_blas_vector_at(y, i) = mc_cadd(mc_blas_vector_at(y, i), mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                   = mc_cadd(temp2, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
				}
				mc_blas_vector_at(y, j) = mc_cadd(mc_blas_vector_at(y, j), mc_cmul(alpha, temp2));
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                    = mc_cmul(alpha, mc_blas_vector_at(x, jx));
				temp2                    = zero;
				mc_blas_vector_at(y, jy) = mc_cadd(mc_blas_vector_at(y, jy), mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, 1, j)));
				l                        = 1 - j;
				ix                       = jx;
				iy                       = jy;
				for (i = (j + 1); i <= mc_minmag(n , j + k); ++i) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_cadd(mc_blas_vector_at(y, iy), mc_cmul(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                    = mc_cadd(temp2, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
				}
				mc_blas_vector_at(y, jy) = mc_cadd(mc_blas_vector_at(y, jy), mc_cmul(alpha, temp2));
				jx                       = jx + incx;
				jy                       = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_qsbmv -

MC_TARGET_FUNC void mc_blas_qsbmv(const char uplo, const int n, const int k, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * a, const int lda, const mc_complex_long_double_t * x, const int incx, const mc_complex_long_double_t beta, mc_complex_long_double_t * y, const int incy)
{
	const mc_complex_long_double_t one = mc_cmplxl(1.0L, 0.0L), zero =  mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (k < 0) {
		info = 3;
	} else if (lda < k + 1) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	} else if (incy == 0) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("QSBMV ", info);
		return;
	}

	if (n == 0 || (mc_ciseql(alpha, zero) && mc_ciseql(beta, one))) {
		return;
	}

	if (incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (n - 1) * incx;
	}
	if (incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (n - 1) * incy;
	}

	if (!mc_ciseql(beta, one)) {
		if (incy == 1) {
			if (mc_ciseql(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = zero;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_cmull(beta, mc_blas_vector_at(y, i));
				}
			}
		} else {
			iy = ky;
			if (mc_ciseql(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = zero;
					iy                       = iy + incy;
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (i = 1; i <= n; ++i) {
					mc_blas_vector_at(y, iy) = mc_cmull(beta, mc_blas_vector_at(y, iy));
					iy                       = iy + incy;
				}
			}
		}
	}

	if (mc_ciseql(alpha, zero)) {
		return;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		kplus1 = k + 1;
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_cmull(alpha, mc_blas_vector_at(x, j));
				temp2 = zero;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_caddl(mc_blas_vector_at(y, i), mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                   = mc_caddl(temp2, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
				}
				mc_blas_vector_at(y, j) = mc_caddl(mc_blas_vector_at(y, j), mc_caddl(mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, kplus1, j)), mc_cmull(alpha, temp2)));
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_cmull(alpha, mc_blas_vector_at(x, jx));
				temp2 = zero;
				ix    = kx;
				iy    = ky;
				l     = kplus1 - j;
				for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
					mc_blas_vector_at(y, iy) = mc_caddl(mc_blas_vector_at(y, iy), mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                    = mc_caddl(temp2, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_caddl(mc_blas_vector_at(y, jy), mc_caddl(mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, kplus1, j)), mc_cmull(alpha, temp2)));
				jx                       = jx + incx;
				jy                       = jy + incy;
				if (j > k) {
					kx = kx + incx;
					ky = ky + incy;
				}
			}
		}
	} else {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                   = mc_cmull(alpha, mc_blas_vector_at(x, j));
				temp2                   = zero;
				mc_blas_vector_at(y, j) = mc_caddl(mc_blas_vector_at(y, j), mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, 1, j)));
				l                       = 1 - j;
				for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
					mc_blas_vector_at(y, i) = mc_caddl(mc_blas_vector_at(y, i), mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                   = mc_caddl(temp2, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
				}
				mc_blas_vector_at(y, j) = mc_caddl(mc_blas_vector_at(y, j), mc_cmull(alpha, temp2));
			}
		} else {
			jx = kx;
			jy = ky;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1                    = mc_cmull(alpha, mc_blas_vector_at(x, jx));
				temp2                    = zero;
				mc_blas_vector_at(y, jy) = mc_caddl(mc_blas_vector_at(y, jy), mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, 1, j)));
				l                        = 1 - j;
				ix                       = jx;
				iy                       = jy;
				for (i = (j + 1); i <= mc_minmag(n , j + k); ++i) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_caddl(mc_blas_vector_at(y, iy), mc_cmull(temp1, mc_blas_matrix_at(a, lda, n, l + i, j)));
					temp2                    = mc_caddl(temp2, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
				}
				mc_blas_vector_at(y, jy) = mc_caddl(mc_blas_vector_at(y, jy), mc_cmull(alpha, temp2));
				jx                       = jx + incx;
				jy                       = jy + incy;
			}
		}
	}
}

#endif /* !MC_BLAS_SBMV_H */

/* EOF */