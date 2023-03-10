//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_syrk.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?syrk performs a rank k operation:
 *    c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c.
 *
 * \synopsis
 *    void ?syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
 *    real-floating alpha, beta
 *    int           k, lda, ldc, n
 *    char          trans, uplo
 *    real-floating a(lda, *), c(ldc, *)
 *
 * \purpose
 *    ?syrk performs a rank k operation: c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c
 *    where alpha and beta are scalars, `c` is an n by n symmetric matrix and `a` is n by k matrix
 *    in the first case and a k by n matrix in the second case.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `c`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `c` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `c` is to be referenced.
 * 
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' c=alpha*a*a' + beta*c + beta*c.
 *    trans='T' or 't' c=alpha*a'*a + beta*c.
 *    trans='C' or 'c' c=alpha*a'*a + beta*c.
 *
 *    [in] n     - int. Specifies the order of the matrix `c`, n must be at least zero.
 *    [in] k     - int. With trans='N' or 'n', k specifies the number of columns of the matrix `a`,
 *    and with trans='T' or 't' or trans='C' or 'c', k specifies the number of rows of the matrix `a`.
 *    k must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, ka), where ka is k when trans='N' or 'n', and is n otherwise.
 *    With trans='T' or 't', the leading n by k part of the array `a` must contain the matrix `a`, otherwise the leading
 *    k by n part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When trans='N' or 'n' then lda must be at least max(1, n),
 *    otherwise lda must be at least max(1, k).
 *
 *    [in] beta  - real-floating. Specifies the scalar beta. When beta is zero then `c` need not be set on input.
 *
 *    [out] c    - real-floating array of dimension (ldc, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `c` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `c` is not referenced. The upper triangular
 *    part of the array `c` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `c` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `c` is not referenced. The lower triangular
 *    part of the array `c` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, n).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Iain S. Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>
#include <macadam/lapack/blas/mc_blas_xerbla.h>
#include <macadam/details/math/mc_cadd.h>
#include <macadam/details/math/mc_ciseq.h>
#include <macadam/details/math/mc_cmul.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_SYRK_H
#define MC_BLAS_SYRK_H

#pragma mark - mc_blas_ssyrk -

MC_TARGET_FUNC void mc_blas_ssyrk(const char uplo, const char trans, const int n, const int k, const float alpha, const float * a, const int lda, const float beta, float * c, const int ldc)
{
	const float one = 1.0f, zero = 0.0f;

	float temp;
	int i, info, j, l, nrowa, ka;
	int upper;

	if (mc_blas_lsame(trans, 'N')) {
		ka    = k;
		nrowa = n;
		mc_unused(ka);
	} else {
		ka    = n;
		nrowa = k;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (n < 0) {
		info = 3;
	} else if (k < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldc < mc_maxmag(1, n)) {
		info = 10;
	}
	if (info != 0) {
		mc_blas_xerbla("SSYRK ", info);
		return;
	}

	if (n == 0 || ((alpha == zero || k == 0) && beta == one)) {
		return;
	}

	if (alpha == zero) {
		if (upper) {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= j; ++i) {
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
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		} else {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = j; i <= n; ++i) {
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
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					if (mc_blas_matrix_at(a, lda, ka, j, l) != zero) {
						temp = alpha * mc_blas_matrix_at(a, lda, ka, j, l);
						for (i = 1; i <= j; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				if (beta == zero) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					if (mc_blas_matrix_at(a, lda, ka, j, l) != zero) {
						temp = alpha * mc_blas_matrix_at(a, lda, ka, j, l);
						for (i = j; i <= n; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
						}
					}
				}
			}
		}
	} else {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= j; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(a, lda, ka, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
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
				for (i = j; i <= n; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(a, lda, ka, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_dsyrk -

MC_TARGET_FUNC void mc_blas_dsyrk(const char uplo, const char trans, const int n, const int k, const double alpha, const double * a, const int lda, const double beta, double * c, const int ldc)
{
	const double one = 1.0, zero = 0.0;

	double temp;
	int i, info, j, l, nrowa, ka;
	int upper;

	if (mc_blas_lsame(trans, 'N')) {
		ka    = k;
		nrowa = n;
		mc_unused(ka);
	} else {
		ka    = n;
		nrowa = k;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (n < 0) {
		info = 3;
	} else if (k < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldc < mc_maxmag(1, n)) {
		info = 10;
	}
	if (info != 0) {
		mc_blas_xerbla("DSYRK ", info);
		return;
	}

	if (n == 0 || ((alpha == zero || k == 0) && beta == one)) {
		return;
	}

	if (alpha == zero) {
		if (upper) {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= j; ++i) {
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
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		} else {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = j; i <= n; ++i) {
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
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					if (mc_blas_matrix_at(a, lda, ka, j, l) != zero) {
						temp = alpha * mc_blas_matrix_at(a, lda, ka, j, l);
						for (i = 1; i <= j; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				if (beta == zero) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					if (mc_blas_matrix_at(a, lda, ka, j, l) != zero) {
						temp = alpha * mc_blas_matrix_at(a, lda, ka, j, l);
						for (i = j; i <= n; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
						}
					}
				}
			}
		}
	} else {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= j; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(a, lda, ka, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
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
				for (i = j; i <= n; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(a, lda, ka, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_lsyrk -

MC_TARGET_FUNC void mc_blas_lsyrk(const char uplo, const char trans, const int n, const int k, const long double alpha, const long double * a, const int lda, const long double beta, long double * c, const int ldc)
{
	const long double one = 1.0L, zero = 0.0L;

	long double temp;
	int i, info, j, l, nrowa, ka;
	int upper;

	if (mc_blas_lsame(trans, 'N')) {
		ka    = k;
		nrowa = n;
		mc_unused(ka);
	} else {
		ka    = n;
		nrowa = k;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (n < 0) {
		info = 3;
	} else if (k < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldc < mc_maxmag(1, n)) {
		info = 10;
	}
	if (info != 0) {
		mc_blas_xerbla("LSYRK ", info);
		return;
	}

	if (n == 0 || ((alpha == zero || k == 0) && beta == one)) {
		return;
	}

	if (alpha == zero) {
		if (upper) {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= j; ++i) {
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
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		} else {
			if (beta == zero) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = j; i <= n; ++i) {
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
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					if (mc_blas_matrix_at(a, lda, ka, j, l) != zero) {
						temp = alpha * mc_blas_matrix_at(a, lda, ka, j, l);
						for (i = 1; i <= j; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				if (beta == zero) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					if (mc_blas_matrix_at(a, lda, ka, j, l) != zero) {
						temp = alpha * mc_blas_matrix_at(a, lda, ka, j, l);
						for (i = j; i <= n; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
						}
					}
				}
			}
		}
	} else {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= j; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(a, lda, ka, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
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
				for (i = j; i <= n; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(a, lda, ka, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
	}
}

/* \name
 *    ?syrk performs a rank k operation:
 *    c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c.
 *
 * \synopsis
 *    void ?syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
 *    complex alpha, beta
 *    int     k, lda, ldc, n
 *    char    trans, uplo
 *    complex a(lda, *), c(ldc, *)
 *
 * \purpose
 *    ?syrk performs a rank k operation: c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c
 *    where alpha and beta are scalars, `c` is an n by n symmetric matrix and `a` is n by k matrix
 *    in the first case and a k by n matrix in the second case.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `c`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `c` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `c` is to be referenced.
 * 
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' c=alpha*a*a' + beta*c + beta*c.
 *    trans='T' or 't' c=alpha*a'*a + beta*c.
 *
 *    [in] n     - int. Specifies the order of the matrix `c`, n must be at least zero.
 *    [in] k     - int. With trans='N' or 'n', k specifies the number of columns of the matrix `a`,
 *    and with trans='T' or 't' or trans='C' or 'c', k specifies the number of rows of the matrix `a`.
 *    k must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, ka), where ka is k when trans='N' or 'n', and is n otherwise.
 *    With trans='T' or 't', the leading n by k part of the array `a` must contain the matrix `a`, otherwise the
 *    leading k by n part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When trans='N' or 'n' then lda must be at least max(1, n),
 *    otherwise lda must be at least max(1, k).
 *
 *    [in] beta  - complex. Specifies the scalar beta. When beta is zero then `c` need not be set on input.
 *
 *    [out] c    - complex array of dimension (ldc, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `c` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `c` is not referenced. The upper triangular
 *    part of the array `c` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `c` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `c` is not referenced. The lower triangular
 *    part of the array `c` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, n).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Iain S. Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#pragma mark - mc_blas_csyrk -

MC_TARGET_FUNC void mc_blas_csyrk(const char uplo, const char trans, const int n, const int k, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t beta, mc_complex_float_t * c, const int ldc)
{
	const mc_complex_float_t one = mc_cmplxf(1.0f, 0.0f), zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp;
	int i, info, j, l, nrowa, ka;
	int upper;

	if (mc_blas_lsame(trans, 'N')) {
		ka    = k;
		nrowa = n;
		mc_unused(ka);
	} else {
		ka    = n;
		nrowa = k;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T')) {
		info = 2;
	} else if (n < 0) {
		info = 3;
	} else if (k < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldc < mc_maxmag(1, n)) {
		info = 10;
	}
	if (info != 0) {
		mc_blas_xerbla("CSYRK ", info);
		return;
	}

	if (n == 0 || ((mc_ciseqf(alpha, zero) || k == 0) && mc_ciseqf(beta, one))) {
		return;
	}

	if (mc_ciseqf(alpha, zero)) {
		if (upper) {
			if (mc_ciseqf(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= j; ++i) {
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
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
			}
		} else {
			if (mc_ciseqf(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = j; i <= n; ++i) {
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
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseqf(beta, zero)) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseqf(beta, one)) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_blas_matrix_at(a, lda, ka, j, l);
					if (!mc_ciseqf(temp, zero)) {
						temp = mc_cmulf(alpha, temp);
						for (i = 1; i <= j; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmulf(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
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
				if (mc_ciseqf(beta, zero)) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseqf(beta, one)) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_blas_matrix_at(a, lda, ka, j, l);
					if (!mc_ciseqf(temp, zero)) {
						temp = mc_cmulf(alpha, temp);
						for (i = j; i <= n; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmulf(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
						}
					}
				}
			}
		}
	} else {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= j; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(a, lda, ka, l, j)));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(alpha, temp), mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
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
				for (i = j; i <= n; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(a, lda, ka, l, j)));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(alpha, temp), mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_zsyrk -

MC_TARGET_FUNC void mc_blas_zsyrk(const char uplo, const char trans, const int n, const int k, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t beta, mc_complex_double_t * c, const int ldc)
{
	const mc_complex_double_t one = mc_cmplx(1.0, 0.0), zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp;
	int i, info, j, l, nrowa, ka;
	int upper;

	if (mc_blas_lsame(trans, 'N')) {
		ka    = k;
		nrowa = n;
		mc_unused(ka);
	} else {
		ka    = n;
		nrowa = k;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T')) {
		info = 2;
	} else if (n < 0) {
		info = 3;
	} else if (k < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldc < mc_maxmag(1, n)) {
		info = 10;
	}
	if (info != 0) {
		mc_blas_xerbla("ZSYRK ", info);
		return;
	}

	if (n == 0 || ((mc_ciseq(alpha, zero) || k == 0) && mc_ciseq(beta, one))) {
		return;
	}

	if (mc_ciseq(alpha, zero)) {
		if (upper) {
			if (mc_ciseq(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= j; ++i) {
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
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
			}
		} else {
			if (mc_ciseq(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = j; i <= n; ++i) {
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
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseq(beta, zero)) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseq(beta, one)) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_blas_matrix_at(a, lda, ka, j, l);
					if (!mc_ciseq(temp, zero)) {
						temp = mc_cmul(alpha, temp);
						for (i = 1; i <= j; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmul(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
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
				if (mc_ciseq(beta, zero)) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseq(beta, one)) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_blas_matrix_at(a, lda, ka, j, l);
					if (!mc_ciseq(temp, zero)) {
						temp = mc_cmul(alpha, temp);
						for (i = j; i <= n; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmul(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
						}
					}
				}
			}
		}
	} else {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= j; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(a, lda, ka, l, j)));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(alpha, temp), mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
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
				for (i = j; i <= n; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(a, lda, ka, l, j)));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(alpha, temp), mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_qsyrk -

MC_TARGET_FUNC void mc_blas_qsyrk(const char uplo, const char trans, const int n, const int k, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * a, const int lda, const mc_complex_long_double_t beta, mc_complex_long_double_t * c, const int ldc)
{
	const mc_complex_long_double_t one = mc_cmplxl(1.0L, 0.0L), zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp;
	int i, info, j, l, nrowa, ka;
	int upper;

	if (mc_blas_lsame(trans, 'N')) {
		ka    = k;
		nrowa = n;
		mc_unused(ka);
	} else {
		ka    = n;
		nrowa = k;
		mc_unused(ka);
	}
	upper = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T')) {
		info = 2;
	} else if (n < 0) {
		info = 3;
	} else if (k < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 7;
	} else if (ldc < mc_maxmag(1, n)) {
		info = 10;
	}
	if (info != 0) {
		mc_blas_xerbla("QSYRK ", info);
		return;
	}

	if (n == 0 || ((mc_ciseql(alpha, zero) || k == 0) && mc_ciseql(beta, one))) {
		return;
	}

	if (mc_ciseql(alpha, zero)) {
		if (upper) {
			if (mc_ciseql(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= j; ++i) {
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
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
			}
		} else {
			if (mc_ciseql(beta, zero)) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = j; i <= n; ++i) {
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
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
			}
		}
		return;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseql(beta, zero)) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseql(beta, one)) {
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_blas_matrix_at(a, lda, ka, j, l);
					if (!mc_ciseql(temp, zero)) {
						temp = mc_cmull(alpha, temp);
						for (i = 1; i <= j; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmull(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
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
				if (mc_ciseql(beta, zero)) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseql(beta, one)) {
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_blas_matrix_at(a, lda, ka, j, l);
					if (!mc_ciseql(temp, zero)) {
						temp = mc_cmull(alpha, temp);
						for (i = j; i <= n; ++i) {
							mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmull(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
						}
					}
				}
			}
		}
	} else {
		if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= j; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(a, lda, ka, l, j)));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(alpha, temp), mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
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
				for (i = j; i <= n; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(a, lda, ka, l, j)));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(alpha, temp), mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	}
}

#endif /* !MC_BLAS_SYRK_H */

/* EOF */