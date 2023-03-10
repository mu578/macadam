//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_gemm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?gemm performs one of the matrix-matrix operations:
 *    c=alpha*op(a)*op(b) + beta*c where op(x)=x or op(x)=x'.
 *
 * \synopsis
 *    void ?gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 *    real-floating alpha, beta
 *    int           k, lda, ldb, ldc, m, n
 *    char          transa, transb
 *    real-floating a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?gemm performs one of the matrix-matrix operations: c=alpha*op(a)*op(b) + beta*c where
 *    op(x)=x or op(x)=x' alpha and beta are scalars, and a, b and c are matrices, with op(a)
 *    an m by k matrix, op(b) a k by n matrix and c an m by n matrix.
 *
 * \parameters
 *    [in] transa - char. Specifies the form of op(a) to be used in the matrix multiplication as follows:
 *    transa='N' or 'n', op(a)=a.
 *    transa='T' or 't', op(a)=a'.
 *    transa='C' or 'c', op(a)=a'.
 *
 *    [in] transb - char. Specifies the form of op(b) to be used in the matrix multiplication as follows:
 *    transb='N' or 'n', op(b)=b.
 *    transb='T' or 't', op(b)=b'.
 *    transb='C' or 'c', op(b)=b'.
 *
 *    [in] m      - int. Specifies the number of rows of the matrix op(a) and of the matrix `c`, m must be
 *    at least zero.
 *
 *    [in] n      - int. Specifies the number of columns of the matrix op(b) and the number of columns of
 *    the matrix `c`, n must be at least zero.
 *
 *    [in] k      - int. Specifies  the number of columns of the matrix op(a) and the number of rows of
 *    the matrix op(b), k must be at least zero.
 *
 *    [in] alpha  - real-floating. Specifies the scalar alpha.
 *
 *    [in] a      - real-floating array of dimension (lda, ka), where ka is k when transa='N' or 'n' and
 *    is m otherwise. Prior entry with transa='N' or 'n', the leading m by k part of the array `a` must
 *    contain the matrix `a`, otherwise the leading k by m part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda    - int. Specifies the first dimension of `a`. When transa='N' or 'n' then
 *    lda must be at least max(1, m), otherwise lda must be at least max(1, k).
 *
 *    [in] b      - real-floating array of dimension (ldb, kb), where kb is n when transb='N' or 'n' and
 *    is k otherwise. Prior entry with transb='N' or 'n', the leading k by n part of the array b must
 *    contain the matrix `b`, otherwise the leading n by k part of the array b must contain the matrix `b`.
 *
 *    [in] ldb    - int. Specifies the first dimension of `b`. When transb='N' or 'n' then
 *    ldb must be at least max(1, k), otherwise ldb must be at least max(1, n).
 *
 *    [in] beta   - real-floating. Specifies the scalar beta. When beta is supplied as zero then c need
 *    not be set on input.
 *
 *    [out] c     - real-floating array of dimension (ldc, n). Prior entry the leading  m by n part of the
 *    array c must contain the matrix `c`, except when beta is set to zero, in which case c need not be set
 *    on entry, c is overwritten by the m by n matrix (alpha*op(a)*op(b) + beta*c).
 *
 *    [in] ldc    - int. Specifies the first dimension of `c`, ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Laboratory.
 *     \author Iain Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>
#include <macadam/lapack/blas/mc_blas_xerbla.h>
#include <macadam/details/math/mc_cadd.h>
#include <macadam/details/math/mc_ciseq.h>
#include <macadam/details/math/mc_conj.h>
#include <macadam/details/math/mc_cmul.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_GEMM_H
#define MC_BLAS_GEMM_H

#pragma mark - mc_blas_sgemm -

MC_TARGET_FUNC void mc_blas_sgemm(const char transa, const char transb, const int m, const int n, const int k, const float alpha, const float * a, const int lda, const float * b, const int ldb, const float beta, float * c, const int ldc)
{
	const float one = 1.0f, zero = 0.0f;

	float temp;
	int i, info, j, l, ncola, nrowa, nrowb, ka, kb;
	int nota, notb;

	nota = mc_blas_lsame(transa, 'N');
	notb = mc_blas_lsame(transb, 'N');

	if (nota) {
		ka    = k;
		nrowa = m;
		ncola = k;
		mc_unused(ka);
		mc_unused(ncola);
	} else {
		ka    = m;
		nrowa = k;
		ncola = m;
		mc_unused(ka);
		mc_unused(ncola);
	}
	if (notb) {
		kb    = n;
		nrowb = k;
		mc_unused(kb);
	} else {
		kb    = k;
		nrowb = n;
		mc_unused(kb);
	}

	info = 0;
	if (!nota && !mc_blas_lsame(transa, 'C') && !mc_blas_lsame(transa, 'T')) {
		info = 1;
	} else if (!notb && !mc_blas_lsame(transb, 'C') && !mc_blas_lsame(transb, 'T')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 8;
	} else if (ldb < mc_maxmag(1, nrowb)) {
		info = 10;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 13;
	}
	if (info != 0) {
		mc_blas_xerbla("SGEMM ", info);
		return;
	}

	if (m == 0 || n == 0 || ((alpha == zero || k == 0) && beta == one)) {
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

	if (notb) {
		if (nota) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = alpha * mc_blas_matrix_at(b, ldb, kb, l, j);
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(b, ldb, kb, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
	} else {
		if (nota) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = alpha * mc_blas_matrix_at(b, ldb, kb, j, l);
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(b, ldb, kb, j, l));
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

#pragma mark - mc_blas_dgemm -

MC_TARGET_FUNC void mc_blas_dgemm(const char transa, const char transb, const int m, const int n, const int k, const double alpha, const double * a, const int lda, const double * b, const int ldb, const double beta, double * c, const int ldc)
{
	const double one = 1.0, zero = 0.0;

	double temp;
	int i, info, j, l, ncola, nrowa, nrowb, ka, kb;
	int nota, notb;

	nota = mc_blas_lsame(transa, 'N');
	notb = mc_blas_lsame(transb, 'N');

	if (nota) {
		ka    = k;
		nrowa = m;
		ncola = k;
		mc_unused(ka);
		mc_unused(ncola);
	} else {
		ka    = m;
		nrowa = k;
		ncola = m;
		mc_unused(ka);
		mc_unused(ncola);
	}
	if (notb) {
		kb    = n;
		nrowb = k;
		mc_unused(kb);
	} else {
		kb    = k;
		nrowb = n;
		mc_unused(kb);
	}

	info = 0;
	if (!nota && !mc_blas_lsame(transa, 'C') && !mc_blas_lsame(transa, 'T')) {
		info = 1;
	} else if (!notb && !mc_blas_lsame(transb, 'C') && !mc_blas_lsame(transb, 'T')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 8;
	} else if (ldb < mc_maxmag(1, nrowb)) {
		info = 10;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 13;
	}
	if (info != 0) {
		mc_blas_xerbla("DGEMM ", info);
		return;
	}

	if (m == 0 || n == 0 || ((alpha == zero || k == 0) && beta == one)) {
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

	if (notb) {
		if (nota) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = alpha * mc_blas_matrix_at(b, ldb, kb, l, j);
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(b, ldb, kb, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
	} else {
		if (nota) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = alpha * mc_blas_matrix_at(b, ldb, kb, j, l);
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(b, ldb, kb, j, l));
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

#pragma mark - mc_blas_lgemm -

MC_TARGET_FUNC void mc_blas_lgemm(const char transa, const char transb, const int m, const int n, const int k, const long double alpha, const long double * a, const int lda, const long double * b, const int ldb, const long double beta, long double * c, const int ldc)
{
	const long double one = 1.0L, zero = 0.0L;

	long double temp;
	int i, info, j, l, ncola, nrowa, nrowb, ka, kb;
	int nota, notb;

	nota = mc_blas_lsame(transa, 'N');
	notb = mc_blas_lsame(transb, 'N');

	if (nota) {
		ka    = k;
		nrowa = m;
		ncola = k;
		mc_unused(ka);
		mc_unused(ncola);
	} else {
		ka    = m;
		nrowa = k;
		ncola = m;
		mc_unused(ka);
		mc_unused(ncola);
	}
	if (notb) {
		kb    = n;
		nrowb = k;
		mc_unused(kb);
	} else {
		kb    = k;
		nrowb = n;
		mc_unused(kb);
	}

	info = 0;
	if (!nota && !mc_blas_lsame(transa, 'C') && !mc_blas_lsame(transa, 'T')) {
		info = 1;
	} else if (!notb && !mc_blas_lsame(transb, 'C') && !mc_blas_lsame(transb, 'T')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 8;
	} else if (ldb < mc_maxmag(1, nrowb)) {
		info = 10;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 13;
	}
	if (info != 0) {
		mc_blas_xerbla("LGEMM ", info);
		return;
	}

	if (m == 0 || n == 0 || ((alpha == zero || k == 0) && beta == one)) {
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

	if (notb) {
		if (nota) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = alpha * mc_blas_matrix_at(b, ldb, kb, l, j);
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(b, ldb, kb, l, j));
					}
					if (beta == zero) {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp;
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = alpha * temp + beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
			}
		}
	} else {
		if (nota) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (beta == zero) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (beta != one) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = beta * mc_blas_matrix_at(c, ldc, n, i, j);
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = alpha * mc_blas_matrix_at(b, ldb, kb, j, l);
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_blas_matrix_at(c, ldc, n, i, j) + (temp * mc_blas_matrix_at(a, lda, ka, i, l));
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
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = temp + (mc_blas_matrix_at(a, lda, ka, l, i) * mc_blas_matrix_at(b, ldb, kb, j, l));
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
 *    ?gemm performs one of the matrix-matrix operations:
 *    c=alpha*op(a)*op(b) + beta*c where op(x)=x or op(x)=x' or op(x)=x_.
 *
 * \synopsis
 *    void ?gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 *    complex alpha, beta
 *    int     k, lda, ldb, ldc, m, n
 *    char    transa, transb
 *    complex a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?gemm performs one of the matrix-matrix operations: c=alpha*op(a)*op(b) + beta*c where
 *    op(x)=x or op(x)=x' or op(x)=x_ alpha and beta are scalars, and a, b and c are matrices,
 *    with op(a) an m by k matrix, op(b) a k by n matrix and c an m by n matrix.
 *
 * \parameters
 *    [in] transa - char. Specifies the form of op(a) to be used in the matrix multiplication as follows:
 *    transa='N' or 'n', op(a)=a.
 *    transa='T' or 't', op(a)=a'.
 *    transa='C' or 'c', op(a)=a_.
 *
 *    [in] transb - char. Specifies the form of op(b) to be used in the matrix multiplication as follows:
 *    transb='N' or 'n', op(b)=b.
 *    transb='T' or 't', op(b)=b'.
 *    transb='C' or 'c', op(b)=b_.
 *
 *    [in] m      - int. Specifies the number of rows of the matrix op(a) and of the matrix `c`, m must be
 *    at least zero.
 *
 *    [in] n      - int. Specifies the number of columns of the matrix op(b) and the number of columns of
 *    the matrix `c`, n must be at least zero.
 *
 *    [in] k      - int. Specifies  the number of columns of the matrix op(a) and the number of rows of
 *    the matrix op(b), k must be at least zero.
 *
 *    [in] alpha  - complex. Specifies the scalar alpha.
 *
 *    [in] a      - complex array of dimension (lda, ka), where ka is k when transa='N' or 'n' and
 *    is m otherwise. Prior entry with transa='N' or 'n', the leading m by k part of the array `a` must
 *    contain the matrix `a`, otherwise the leading k by m part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda    - int. Specifies the first dimension of `a`. When transa='N' or 'n' then
 *    lda must be at least max(1, m), otherwise lda must be at least max(1, k).
 *
 *    [in] b      - complex array of dimension (ldb, kb), where kb is n when transb='N' or 'n' and
 *    is k otherwise. Prior entry with transb='N' or 'n', the leading k by n part of the array b must
 *    contain the matrix `b`, otherwise the leading n by k part of the array b must contain the matrix `b`.
 *
 *    [in] ldb    - int. Specifies the first dimension of `b`. When transb='N' or 'n' then
 *    ldb must be at least max(1, k), otherwise ldb must be at least max(1, n).
 *
 *    [in] beta   - complex. Specifies the scalar beta. When beta is supplied as zero then c need
 *    not be set on input.
 *
 *    [out] c     - complex array of dimension (ldc, n). Prior entry the leading  m by n part of the
 *    array c must contain the matrix `c`, except when beta is set to zero, in which case c need not be set
 *    on entry, c is overwritten by the m by n matrix (alpha*op(a)*op(b) + beta*c).
 *
 *    [in] ldc    - int. Specifies the first dimension of `c`, ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Laboratory.
 *     \author Iain Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#pragma mark - mc_blas_cgemm -

MC_TARGET_FUNC void mc_blas_cgemm(const char transa, const char transb, const int m, const int n, const int k, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * b, const int ldb, const mc_complex_float_t beta, mc_complex_float_t * c, const int ldc)
{
	const mc_complex_float_t one = mc_cmplxf(1.0f, 0.0f), zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp;
	int i, info, j, l, ncola, nrowa, nrowb, ka, kb;
	int conja, conjb, nota, notb;

	nota  = mc_blas_lsame(transa, 'N');
	notb  = mc_blas_lsame(transb, 'N');
	conja = mc_blas_lsame(transa, 'C');
	conjb = mc_blas_lsame(transb, 'C');

	if (nota) {
		ka    = k;
		nrowa = m;
		ncola = k;
		mc_unused(ka);
		mc_unused(ncola);
	} else {
		ka    = m;
		nrowa = k;
		ncola = m;
		mc_unused(ka);
		mc_unused(ncola);
	}
	if (notb) {
		kb    = n;
		nrowb = k;
		mc_unused(kb);
	} else {
		kb    = k;
		nrowb = n;
		mc_unused(kb);
	}

	info = 0;
	if (!nota && !conja && !mc_blas_lsame(transa, 'T')) {
		info = 1;
	} else if (!notb && !conjb && !mc_blas_lsame(transb, 'T')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 8;
	} else if (ldb < mc_maxmag(1, nrowb)) {
		info = 10;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 13;
	}
	if (info != 0) {
		mc_blas_xerbla("CGEMM ", info);
		return;
	}

	if (m == 0 || n == 0 || ((mc_ciseqf(alpha, zero) || k == 0) && mc_ciseqf(beta, one))) {
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

	if (notb) {
		if (nota) {
//!# c=alpha*a*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseqf(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseqf(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmulf(alpha, mc_blas_matrix_at(b, ldb, kb, l, j));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmulf(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		} else if (conja) {
//!# c=alpha*a_*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_matrix_at(a, lda, ka, l, i)), mc_blas_matrix_at(b, ldb, kb, l, j)));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(alpha, temp), mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a'*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(b, ldb, kb, l, j)));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(alpha, temp), mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	} else if (nota) {
		if (conjb) {
//!# c=alpha*a*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseqf(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseqf(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmulf(alpha, mc_conjf(mc_blas_matrix_at(b, ldb, kb, j, l)));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmulf(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		} else {
//!# c=alpha*a*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseqf(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseqf(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmulf(alpha, mc_blas_matrix_at(b, ldb, kb, j, l));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmulf(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		}
	} else if (conja) {
		if (conjb) {
//!# c=alpha*a_*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_matrix_at(a, lda, ka, l, i)), mc_conjf(mc_blas_matrix_at(b, ldb, kb, j, l))));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(alpha, temp), mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a_*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_matrix_at(a, lda, ka, l, i)), mc_blas_matrix_at(b, ldb, kb, j, l)));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(alpha, temp), mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	} else {
		if (conjb) {
//!# c=alpha*a'*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, ka, l, i), mc_conjf(mc_blas_matrix_at(b, ldb, kb, j, l))));
					}
					if (mc_ciseqf(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmulf(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddf(mc_cmulf(alpha, temp), mc_cmulf(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a'*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(b, ldb, kb, j, l)));
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

#pragma mark - mc_blas_zgemm -

MC_TARGET_FUNC void mc_blas_zgemm(const char transa, const char transb, const int m, const int n, const int k, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * b, const int ldb, const mc_complex_double_t beta, mc_complex_double_t * c, const int ldc)
{
	const mc_complex_double_t one = mc_cmplx(1.0, 0.0), zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp;
	int i, info, j, l, ncola, nrowa, nrowb, ka, kb;
	int conja, conjb, nota, notb;

	nota  = mc_blas_lsame(transa, 'N');
	notb  = mc_blas_lsame(transb, 'N');
	conja = mc_blas_lsame(transa, 'C');
	conjb = mc_blas_lsame(transb, 'C');

	if (nota) {
		ka    = k;
		nrowa = m;
		ncola = k;
		mc_unused(ka);
		mc_unused(ncola);
	} else {
		ka    = m;
		nrowa = k;
		ncola = m;
		mc_unused(ka);
		mc_unused(ncola);
	}
	if (notb) {
		kb    = n;
		nrowb = k;
		mc_unused(kb);
	} else {
		kb    = k;
		nrowb = n;
		mc_unused(kb);
	}

	info = 0;
	if (!nota && !conja && !mc_blas_lsame(transa, 'T')) {
		info = 1;
	} else if (!notb && !conjb && !mc_blas_lsame(transb, 'T')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 8;
	} else if (ldb < mc_maxmag(1, nrowb)) {
		info = 10;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 13;
	}
	if (info != 0) {
		mc_blas_xerbla("ZGEMM ", info);
		return;
	}

	if (m == 0 || n == 0 || ((mc_ciseq(alpha, zero) || k == 0) && mc_ciseq(beta, one))) {
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

	if (notb) {
		if (nota) {
//!# c=alpha*a*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseq(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseq(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmul(alpha, mc_blas_matrix_at(b, ldb, kb, l, j));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmul(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		} else if (conja) {
//!# c=alpha*a_*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_matrix_at(a, lda, ka, l, i)), mc_blas_matrix_at(b, ldb, kb, l, j)));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(alpha, temp), mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a'*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(b, ldb, kb, l, j)));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(alpha, temp), mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	} else if (nota) {
		if (conjb) {
//!# c=alpha*a*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseq(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseq(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmul(alpha, mc_conj(mc_blas_matrix_at(b, ldb, kb, j, l)));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmul(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		} else {
//!# c=alpha*a*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseq(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseq(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmul(alpha, mc_blas_matrix_at(b, ldb, kb, j, l));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmul(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		}
	} else if (conja) {
		if (conjb) {
//!# c=alpha*a_*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_matrix_at(a, lda, ka, l, i)), mc_conj(mc_blas_matrix_at(b, ldb, kb, j, l))));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(alpha, temp), mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a_*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_matrix_at(a, lda, ka, l, i)), mc_blas_matrix_at(b, ldb, kb, j, l)));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(alpha, temp), mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	} else {
		if (conjb) {
//!# c=alpha*a'*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, ka, l, i), mc_conj(mc_blas_matrix_at(b, ldb, kb, j, l))));
					}
					if (mc_ciseq(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmul(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cadd(mc_cmul(alpha, temp), mc_cmul(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a'*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(b, ldb, kb, j, l)));
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

#pragma mark - mc_blas_qgemm -

MC_TARGET_FUNC void mc_blas_qgemm(const char transa, const char transb, const int m, const int n, const int k, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * a, const int lda, const mc_complex_long_double_t * b, const int ldb, const mc_complex_long_double_t beta, mc_complex_long_double_t * c, const int ldc)
{
	const mc_complex_long_double_t one = mc_cmplxl(1.0L, 0.0L), zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp;
	int i, info, j, l, ncola, nrowa, nrowb, ka, kb;
	int conja, conjb, nota, notb;

	nota  = mc_blas_lsame(transa, 'N');
	notb  = mc_blas_lsame(transb, 'N');
	conja = mc_blas_lsame(transa, 'C');
	conjb = mc_blas_lsame(transb, 'C');

	if (nota) {
		ka    = k;
		nrowa = m;
		ncola = k;
		mc_unused(ka);
		mc_unused(ncola);
	} else {
		ka    = m;
		nrowa = k;
		ncola = m;
		mc_unused(ka);
		mc_unused(ncola);
	}
	if (notb) {
		kb    = n;
		nrowb = k;
		mc_unused(kb);
	} else {
		kb    = k;
		nrowb = n;
		mc_unused(kb);
	}

	info = 0;
	if (!nota && !conja && !mc_blas_lsame(transa, 'T')) {
		info = 1;
	} else if (!notb && !conjb && !mc_blas_lsame(transb, 'T')) {
		info = 2;
	} else if (m < 0) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 8;
	} else if (ldb < mc_maxmag(1, nrowb)) {
		info = 10;
	} else if (ldc < mc_maxmag(1, m)) {
		info = 13;
	}
	if (info != 0) {
		mc_blas_xerbla("QGEMM ", info);
		return;
	}

	if (m == 0 || n == 0 || ((mc_ciseql(alpha, zero) || k == 0) && mc_ciseql(beta, one))) {
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

	if (notb) {
		if (nota) {
//!# c=alpha*a*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseql(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseql(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmull(alpha, mc_blas_matrix_at(b, ldb, kb, l, j));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmull(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		} else if (conja) {
//!# c=alpha*a_*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_matrix_at(a, lda, ka, l, i)), mc_blas_matrix_at(b, ldb, kb, l, j)));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(alpha, temp), mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a'*b+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(b, ldb, kb, l, j)));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(alpha, temp), mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	} else if (nota) {
		if (conjb) {
//!# c=alpha*a*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseql(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseql(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmull(alpha, mc_conjl(mc_blas_matrix_at(b, ldb, kb, j, l)));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmull(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		} else {
//!# c=alpha*a*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_ciseql(beta, zero)) {
					for (i = 1; i <= m; ++i) {
						 mc_blas_matrix_at(c, ldc, n, i, j) = zero;
					}
				} else if (!mc_ciseql(beta, one)) {
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j));
					}
				}
				for (l = 1; l <= k; ++l) {
					temp = mc_cmull(alpha, mc_blas_matrix_at(b, ldb, kb, j, l));
					for (i = 1; i <= m; ++i) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_blas_matrix_at(c, ldc, n, i, j), mc_cmull(temp, mc_blas_matrix_at(a, lda, ka, i, l)));
					}
				}
			}
		}
	} else if (conja) {
		if (conjb) {
//!# c=alpha*a_*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_matrix_at(a, lda, ka, l, i)), mc_conjl(mc_blas_matrix_at(b, ldb, kb, j, l))));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(alpha, temp), mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a_*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_matrix_at(a, lda, ka, l, i)), mc_blas_matrix_at(b, ldb, kb, j, l)));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(alpha, temp), mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		}
	} else {
		if (conjb) {
//!# c=alpha*a'*b_+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, ka, l, i), mc_conjl(mc_blas_matrix_at(b, ldb, kb, j, l))));
					}
					if (mc_ciseql(beta, zero)) {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_cmull(alpha, temp);
					} else {
						mc_blas_matrix_at(c, ldc, n, i, j) = mc_caddl(mc_cmull(alpha, temp), mc_cmull(beta, mc_blas_matrix_at(c, ldc, n, i, j)));
					}
				}
			}
		} else {
//!# c=alpha*a'*b'+beta*c
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				for (i = 1; i <= m; ++i) {
					temp = zero;
					for (l = 1; l <= k; ++l) {
						temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, ka, l, i), mc_blas_matrix_at(b, ldb, kb, j, l)));
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

#endif /* !MC_BLAS_GEMM_H */

/* EOF */