//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_symm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?symm performs the matrix-matrix operation:
 *    c=alpha*a*b + beta*c or c=alpha*b*a + beta*c.
 *
 * \synopsis
 *    void ?symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 *    real-floating alpha, beta
 *    int           lda, ldb, ldc, m, n
 *    char          side, uplo
 *    real-floating a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?symm performs the matrix-matrix operation: c=alpha*a*b + beta*c or or c=alpha*b*a + beta*c
 *    where alpha and beta are scalars, `a` is a symmetric matrix and `b` and `c` are m by n matrices.
 *
 * \parameters
 *    [in] side  - char. Specifies whether  the symmetric matrix `a` appears on the left or right
 *    in the  operation as follows:
 *    side='L' or 'l' c=alpha*a*b + beta*c.
 *    side='R' or 'r' c=alpha*a*b + beta*c.
 *
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the symmetric
 *    matrix `a` is to be referenced as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] m     - int. Specifies the number of rows of the matrix `c`, m must be at least zero.
 *    [in] n     - int. Specifies the number of columns of the matrix `c`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, ka), where ka is m when side='L' or 'l' and is n otherwise.
 *    The m by m part of the array `a` must contain the symmetric matrix, such that when uplo='U' or 'u', the leading
 *    m by m upper triangular part of the array `a` must contain the upper triangular part of the symmetric matrix and
 *    the strictly lower triangular part of `a` is not referenced, and when uplo='L' or 'l', the leading m by m lower
 *    triangular part of the array `a` must contain the lower triangular part of the symmetric matrix and the strictly
 *    upper triangular part of `a` is not referenced.
 *
 *    With side='R' or 'r', the n by n part of the array `a` must contain the symmetric matrix, such that when
 *    uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced, and when
 *    uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When side='L' or 'l' then lda must be at least max(1, m),
 *    otherwise lda must be at least max(1, n).
 *
 *    [int] b    - real-floating array of dimension (ldb, n).
 *    The leading m by n part of the array `b` must contain the matrix `b`.
 *
 *    [in] ldb   - int. Specifies the first dimension of `b`. ldb must be at least max(1, m).
 *
 *    [in] beta  - real-floating. Specifies the scalar beta. When beta is zero then `c` need not be set on input.

 *    [out] c    - real-floating array of dimension (ldc, n).
 *    The leading m by n part of the array `c` must contain the matrix `c`, except when beta is zero. The array `c` is
 *    overwritten by the m by n updated resulted matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
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

#ifndef MC_BLAS_SYMM_H
#define MC_BLAS_SYMM_H

#pragma mark - mc_blas_ssymm -

MC_TARGET_FUNC void mc_blas_ssymm(const char side, const char uplo, const int m, const int n, const float alpha, const float * a, const int lda, const float * b, const int ldb, const float beta, float * c, const int ldc)
{
	const float one = 1.0f, zero = 0.0f;

	float temp1, temp2;
	int i, info, j, k, nrowa, ka;
	int upper;

	if (mc_blas_lsame(side, 'L')) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!mc_blas_lsame(side, 'L') && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 9;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 12;
	}
	if (info != 0) {
		mc_blas_xerbla("SSYMM ", info);
		return;
	}

	if (m == 0 || n == 0 || (alpha == zero && beta == one)) {
		return;
	}

	if (alpha == zero) {
		if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = zero;
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp1 = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
					temp2 = zero;
					for (k = 1; k <= (i - 1); ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_blas_matrix_at(c, ldc, n, k, j) + (temp1 * mc_blas_matrix_at(a, lda, ka, k, i));
						temp2                              = temp2 + (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, k, i));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					}
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = m; i >= 1; --i) {
					temp1 = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
					temp2 = zero;
					for (k = i + 1; k <= m; ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_blas_matrix_at(c, ldc, n, k, j) + (temp1 * mc_blas_matrix_at(a, lda, ka, k, i));
						temp2                              = temp2 + (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, k, i));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					}
				}
			}
		}
	} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, j);
			if (beta == zero) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(b, ldb, n, i, j);
				}
			} else {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(b, ldb, n, i, j);
				}
			}
			for (k = 1; k <= (j - 1); ++k) {
				if (upper) {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, k, j);
				} else {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, k);
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp1 * mc_blas_matrix_at(b, ldb, n, i, k));
				}
			}
			for (k = j + 1; k <= n; ++k) {
				if (upper) {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, k);
				} else {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, k, j);
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp1 * mc_blas_matrix_at(b, ldb, n, i, k));
				}
			}
		}
	}
}

#pragma mark - mc_blas_dsymm -

MC_TARGET_FUNC void mc_blas_dsymm(const char side, const char uplo, const int m, const int n, const double alpha, const double * a, const int lda, const double * b, const int ldb, const double beta, double * c, const int ldc)
{
	const double one = 1.0, zero = 0.0;

	double temp1, temp2;
	int i, info, j, k, nrowa, ka;
	int upper;

	if (mc_blas_lsame(side, 'L')) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!mc_blas_lsame(side, 'L') && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 9;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 12;
	}
	if (info != 0) {
		mc_blas_xerbla("DSYMM ", info);
		return;
	}

	if (m == 0 || n == 0 || (alpha == zero && beta == one)) {
		return;
	}

	if (alpha == zero) {
		if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = zero;
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp1 = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
					temp2 = zero;
					for (k = 1; k <= (i - 1); ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_blas_matrix_at(c, ldc, n, k, j) + (temp1 * mc_blas_matrix_at(a, lda, ka, k, i));
						temp2                              = temp2 + (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, k, i));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					}
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = m; i >= 1; --i) {
					temp1 = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
					temp2 = zero;
					for (k = i + 1; k <= m; ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_blas_matrix_at(c, ldc, n, k, j) + (temp1 * mc_blas_matrix_at(a, lda, ka, k, i));
						temp2                              = temp2 + (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, k, i));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					}
				}
			}
		}
	} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, j);
			if (beta == zero) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(b, ldb, n, i, j);
				}
			} else {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(b, ldb, n, i, j);
				}
			}
			for (k = 1; k <= (j - 1); ++k) {
				if (upper) {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, k, j);
				} else {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, k);
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp1 * mc_blas_matrix_at(b, ldb, n, i, k));
				}
			}
			for (k = j + 1; k <= n; ++k) {
				if (upper) {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, k);
				} else {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, k, j);
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp1 * mc_blas_matrix_at(b, ldb, n, i, k));
				}
			}
		}
	}
}

#pragma mark - mc_blas_lsymm -

MC_TARGET_FUNC void mc_blas_lsymm(const char side, const char uplo, const int m, const int n, const long double alpha, const long double * a, const int lda, const long double * b, const int ldb, const long double beta, long double * c, const int ldc)
{
	const long double one = 1.0L, zero = 0.0L;

	long double temp1, temp2;
	int i, info, j, k, nrowa, ka;
	int upper;

	if (mc_blas_lsame(side, 'L')) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!mc_blas_lsame(side, 'L') && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 9;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 12;
	}
	if (info != 0) {
		mc_blas_xerbla("LSYMM ", info);
		return;
	}

	if (m == 0 || n == 0 || (alpha == zero && beta == one)) {
		return;
	}

	if (alpha == zero) {
		if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = zero;
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp1 = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
					temp2 = zero;
					for (k = 1; k <= (i - 1); ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_blas_matrix_at(c, ldc, n, k, j) + (temp1 * mc_blas_matrix_at(a, lda, ka, k, i));
						temp2                              = temp2 + (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, k, i));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					}
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = m; i >= 1; --i) {
					temp1 = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
					temp2 = zero;
					for (k = i + 1; k <= m; ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_blas_matrix_at(c, ldc, n, k, j) + (temp1 * mc_blas_matrix_at(a, lda, ka, k, i));
						temp2                              = temp2 + (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, k, i));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(a, lda, ka, i, i) + alpha * temp2;
					}
				}
			}
		}
	} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, j);
			if (beta == zero) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = temp1 * mc_blas_matrix_at(b, ldb, n, i, j);
				}
			} else {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j) + temp1 * mc_blas_matrix_at(b, ldb, n, i, j);
				}
			}
			for (k = 1; k <= (j - 1); ++k) {
				if (upper) {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, k, j);
				} else {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, k);
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp1 * mc_blas_matrix_at(b, ldb, n, i, k));
				}
			}
			for (k = j + 1; k <= n; ++k) {
				if (upper) {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, j, k);
				} else {
					temp1 = alpha * mc_blas_matrix_at(a, lda, ka, k, j);
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp1 * mc_blas_matrix_at(b, ldb, n, i, k));
				}
			}
		}
	}
}

/* \name
 *    ?symm performs the matrix-matrix operation:
 *    c=alpha*a*b + beta*c or c=alpha*b*a + beta*c.
 *
 * \synopsis
 *    void ?symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 *    complex alpha, beta
 *    int     lda, ldb, ldc, m, n
 *    char    side, uplo
 *    complex a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?symm performs the matrix-matrix operation: c=alpha*a*b + beta*c or c=alpha*b*a + beta*c where
 *    alpha and beta are scalars, `a` is a symmetric matrix and `b` and `c` are m by n matrices.
 *
 * \parameters
 *    [in] side  - char. Specifies whether  the symmetric matrix `a` appears on the left or right
 *    in the  operation as follows:
 *    side='L' or 'l' c=alpha*a*b + beta*c.
 *    side='R' or 'r' c=alpha*a*b + beta*c.
 *
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the symmetric
 *    matrix `a` is to be referenced as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] m     - int. Specifies the number of rows of the matrix `c`, m must be at least zero.
 *    [in] n     - int. Specifies the number of columns of the matrix `c`, n must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, ka), where ka is m when side='L' or 'l' and is n otherwise.
 *    The m by m part of the array `a` must contain the symmetric matrix, such that when uplo='U' or 'u', the leading
 *    m by m upper triangular part of the array `a` must contain the upper triangular part of the symmetric matrix and
 *    the strictly lower triangular part of `a` is not referenced, and when uplo='L' or 'l', the leading m by m lower
 *    triangular part of the array `a` must contain the lower triangular part of the symmetric matrix and the strictly
 *    upper triangular part of `a` is not referenced.
 *
 *    With side='R' or 'r', the n by n part of the array `a` must contain the symmetric matrix, such that when
 *    uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced, and when
 *    uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When side='L' or 'l' then lda must be at least max(1, m),
 *    otherwise lda must be at least max(1, n).
 *
 *    [int] b    - complex array of dimension (lda, n).
 *    The leading m by n part of the array `b` must contain the matrix `b`.
 *
 *    [in] ldb   - int. Specifies the first dimension of `b`. ldb must be at least max(1, m).
 *
 *    [in] beta  - complex. Specifies the scalar beta. When beta is zero then `c` need not be set on input.

 *    [out] c    - complex array of dimension (ldc, n).
 *    The leading m by n part of the array `c` must contain the matrix `c`, except when beta is zero. The array `c` is
 *    overwritten by the m by n updated resulted matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Jeremy Du Croz, Nag Central Office.
 *     \author Sven Hammarling, Nag Central Office.
 *     \author Richard Hanson, Sandia National Labs.
 */

#pragma mark - mc_blas_csymm -

MC_TARGET_FUNC void mc_blas_csymm(const char side, const char uplo, const int m, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * b, const int ldb, const mc_complex_float_t beta, mc_complex_float_t * c, const int ldc)
{
	const mc_complex_float_t one = mc_cmplxf(1.0f, 0.0f), zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp1, temp2;
	int i, info, j, k, nrowa, ka;
	int upper;

	if (mc_blas_lsame(side, 'L')) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!mc_blas_lsame(side, 'L') && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 9;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 12;
	}
	if (info != 0) {
		mc_blas_xerbla("CSYMM ", info);
		return;
	}

	if (m == 0 || n == 0 || (mc_ciseqf(alpha, zero) && mc_ciseqf(beta, one))) {
		return;
	}

	if (mc_ciseqf(alpha, zero)) {
		if (mc_ciseqf(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = zero;
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp1 = mc_cmulf(alpha, mc_blas_matrix_at(b, ldb, n, i, j));
					temp2 = zero;
					for (k = 1; k <= (i - 1); ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, k, j), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, ka, k, i)));
						temp2                              = mc_caddf(temp2, mc_cmulf(mc_blas_matrix_at(b, ldb, n, k, j), mc_blas_matrix_at(a, lda, ka, k, i)));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmulf(alpha, temp2));
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_caddf(mc_cmulf(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmulf(alpha, temp2)));
					}
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = m; i >= 1; --i) {
					temp1 = mc_cmulf(alpha, mc_blas_matrix_at(b, ldb, n, i, j));
					temp2 = zero;
					for (k = i + 1; k <= m; ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, k, j), mc_cmulf(temp1, mc_blas_matrix_at(a, lda, ka, k, i)));
						temp2                              = mc_caddf(temp2, mc_cmulf(mc_blas_matrix_at(b, ldb, n, k, j), mc_blas_matrix_at(a, lda, ka, k, i)));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmulf(alpha, temp2));
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_caddf(mc_cmulf(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmulf(alpha, temp2)));
					}
				}
			}
		}
	} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			temp1 = mc_cmulf(alpha, mc_blas_matrix_at(a, lda, ka, j, j));
			if (mc_ciseqf(beta, zero)) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(temp1, mc_blas_matrix_at(b, ldb, n, i, j));
				}
			} else {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_cmulf(temp1, mc_blas_matrix_at(b, ldb, n, i, j)));
				}
			}
			for (k = 1; k <= (j - 1); ++k) {
				if (upper) {
					temp1 = mc_cmulf(alpha, mc_blas_matrix_at(a, lda, ka, k, j));
				} else {
					temp1 = mc_cmulf(alpha, mc_blas_matrix_at(a, lda, ka, j, k));
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmulf(temp1, mc_blas_matrix_at(b, ldb, n, i, k)));
				}
			}
			for (k = j + 1; k <= n; ++k) {
				if (upper) {
					temp1 = mc_cmulf(alpha, mc_blas_matrix_at(a, lda, ka, j, k));
				} else {
					temp1 = mc_cmulf(alpha, mc_blas_matrix_at(a, lda, ka, k, j));
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmulf(temp1, mc_blas_matrix_at(b, ldb, n, i, k)));
				}
			}
		}
	}
}

#pragma mark - mc_blas_zsymm -

MC_TARGET_FUNC void mc_blas_zsymm(const char side, const char uplo, const int m, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * b, const int ldb, const mc_complex_double_t beta, mc_complex_double_t * c, const int ldc)
{
	const mc_complex_double_t one = mc_cmplx(1.0, 0.0), zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp1, temp2;
	int i, info, j, k, nrowa, ka;
	int upper;

	if (mc_blas_lsame(side, 'L')) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!mc_blas_lsame(side, 'L') && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 9;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 12;
	}
	if (info != 0) {
		mc_blas_xerbla("ZSYMM ", info);
		return;
	}

	if (m == 0 || n == 0 || (mc_ciseq(alpha, zero) && mc_ciseq(beta, one))) {
		return;
	}

	if (mc_ciseq(alpha, zero)) {
		if (mc_ciseq(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = zero;
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp1 = mc_cmul(alpha, mc_blas_matrix_at(b, ldb, n, i, j));
					temp2 = zero;
					for (k = 1; k <= (i - 1); ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, k, j), mc_cmul(temp1, mc_blas_matrix_at(a, lda, ka, k, i)));
						temp2                              = mc_cadd(temp2, mc_cmul(mc_blas_matrix_at(b, ldb, n, k, j), mc_blas_matrix_at(a, lda, ka, k, i)));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmul(alpha, temp2));
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_cadd(mc_cmul(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmul(alpha, temp2)));
					}
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = m; i >= 1; --i) {
					temp1 = mc_cmul(alpha, mc_blas_matrix_at(b, ldb, n, i, j));
					temp2 = zero;
					for (k = i + 1; k <= m; ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, k, j), mc_cmul(temp1, mc_blas_matrix_at(a, lda, ka, k, i)));
						temp2                              = mc_cadd(temp2, mc_cmul(mc_blas_matrix_at(b, ldb, n, k, j), mc_blas_matrix_at(a, lda, ka, k, i)));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmul(alpha, temp2));
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_cadd(mc_cmul(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmul(alpha, temp2)));
					}
				}
			}
		}
	} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			temp1 = mc_cmul(alpha, mc_blas_matrix_at(a, lda, ka, j, j));
			if (mc_ciseq(beta, zero)) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(temp1, mc_blas_matrix_at(b, ldb, n, i, j));
				}
			} else {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_cmul(temp1, mc_blas_matrix_at(b, ldb, n, i, j)));
				}
			}
			for (k = 1; k <= (j - 1); ++k) {
				if (upper) {
					temp1 = mc_cmul(alpha, mc_blas_matrix_at(a, lda, ka, k, j));
				} else {
					temp1 = mc_cmul(alpha, mc_blas_matrix_at(a, lda, ka, j, k));
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmul(temp1, mc_blas_matrix_at(b, ldb, n, i, k)));
				}
			}
			for (k = j + 1; k <= n; ++k) {
				if (upper) {
					temp1 = mc_cmul(alpha, mc_blas_matrix_at(a, lda, ka, j, k));
				} else {
					temp1 = mc_cmul(alpha, mc_blas_matrix_at(a, lda, ka, k, j));
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmul(temp1, mc_blas_matrix_at(b, ldb, n, i, k)));
				}
			}
		}
	}
}

#pragma mark - mc_blas_qsymm -

MC_TARGET_FUNC void mc_blas_qsymm(const char side, const char uplo, const int m, const int n, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * a, const int lda, const mc_complex_long_double_t * b, const int ldb, const mc_complex_long_double_t beta, mc_complex_long_double_t * c, const int ldc)
{
	const mc_complex_long_double_t one = mc_cmplxl(1.0L, 0.0L), zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp1, temp2;
	int i, info, j, k, nrowa, ka;
	int upper;

	if (mc_blas_lsame(side, 'L')) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!mc_blas_lsame(side, 'L') && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 9;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 12;
	}
	if (info != 0) {
		mc_blas_xerbla("QSYMM ", info);
		return;
	}

	if (m == 0 || n == 0 || (mc_ciseql(alpha, zero) && mc_ciseql(beta, one))) {
		return;
	}

	if (mc_ciseql(alpha, zero)) {
		if (mc_ciseql(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = zero;
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp1 = mc_cmull(alpha, mc_blas_matrix_at(b, ldb, n, i, j));
					temp2 = zero;
					for (k = 1; k <= (i - 1); ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, k, j), mc_cmull(temp1, mc_blas_matrix_at(a, lda, ka, k, i)));
						temp2                              = mc_caddl(temp2, mc_cmull(mc_blas_matrix_at(b, ldb, n, k, j), mc_blas_matrix_at(a, lda, ka, k, i)));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmull(alpha, temp2));
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_caddl(mc_cmull(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmull(alpha, temp2)));
					}
				}
			}
		} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = m; i >= 1; --i) {
					temp1 = mc_cmull(alpha, mc_blas_matrix_at(b, ldb, n, i, j));
					temp2 = zero;
					for (k = i + 1; k <= m; ++k) {
						mc_blas_matrix_at(c, ldc, n, k, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, k, j), mc_cmull(temp1, mc_blas_matrix_at(a, lda, ka, k, i)));
						temp2                              = mc_caddl(temp2, mc_cmull(mc_blas_matrix_at(b, ldb, n, k, j), mc_blas_matrix_at(a, lda, ka, k, i)));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmull(alpha, temp2));
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_caddl(mc_cmull(temp1, mc_blas_matrix_at(a, lda, ka, i, i)), mc_cmull(alpha, temp2)));
					}
				}
			}
		}
	} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (j = 1; j <= n; ++j) {
			temp1 = mc_cmull(alpha, mc_blas_matrix_at(a, lda, ka, j, j));
			if (mc_ciseql(beta, zero)) {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(temp1, mc_blas_matrix_at(b, ldb, n, i, j));
				}
			} else {
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)), mc_cmull(temp1, mc_blas_matrix_at(b, ldb, n, i, j)));
				}
			}
			for (k = 1; k <= (j - 1); ++k) {
				if (upper) {
					temp1 = mc_cmull(alpha, mc_blas_matrix_at(a, lda, ka, k, j));
				} else {
					temp1 = mc_cmull(alpha, mc_blas_matrix_at(a, lda, ka, j, k));
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmull(temp1, mc_blas_matrix_at(b, ldb, n, i, k)));
				}
			}
			for (k = j + 1; k <= n; ++k) {
				if (upper) {
					temp1 = mc_cmull(alpha, mc_blas_matrix_at(a, lda, ka, j, k));
				} else {
					temp1 = mc_cmull(alpha, mc_blas_matrix_at(a, lda, ka, k, j));
				}
				for (i = 1; i <= m; ++i) {
					mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmull(temp1, mc_blas_matrix_at(b, ldb, n, i, k)));
				}
			}
		}
	}
}

#endif /* !MC_BLAS_SYMM_H */

/* EOF */