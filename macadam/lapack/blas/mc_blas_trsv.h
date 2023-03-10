//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_trsv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?trsv Solves a system of linear equations whose coefficients are in a triangular matrix.
 *
 * \synopsis
 *    void ?trsv(uplo, trans, diag, n, a, lda, x, incx)
 *    real-floating alpha, beta
 *    int           incx, lda, n
 *    char          uplo, trans, diag
 *    real-floating a(lda,*), x(*)
 *
 * \purpose
 *    ?trsv solves one of the systems of equations: x*x=b or a'*x=b where b and x are n element
 *    vectors and `a` is an n by n unit or non-unit, upper or lower triangular matrix. However,
 *    no test for singularity or near-singularity is included in this routine. Such tests must
 *    be performed before calling this routine.
 *
 *    The term b used in the systems of equations listed above represents the right-hand side of
 *    the system. It is important to note that in these subroutines the right-hand side of the
 *    equation is actually provided in the input-output argument x.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the matrix is an upper or lower triangular matrix as follows:
 *    uplo='U' or 'u' a is an upper triangular matrix.
 *    uplo='L' or 'l' a is a lower triangular matrix.
 *
 *    [in] trans - char. Specifies the equations to be solved as follows:
 *    trans='N' or 'n' a*x  = b.
 *    trans='T' or 't' a'*x = b.
 *    trans='C' or 'c' a'*x = b.
 *
 *    [in] diag  - char. Specifies whether or not `a` is unit triangular as follows:
 *    diag='U' or 'u' `a` is assumed to be unit triangular.
 *    diag='N' or 'n' `a` is not assumed to be unit triangular.
 *
 *    [in] n     - int. Specifies the order of the matrix `a`, n must be at least zero; n>=0 and n<=lda
 *
 *    [in] a     - real-floating array of dimension (lda, n). With uplo='U' or 'u', the leading n by n
 *    upper triangular part of the array `a` must contain the upper triangular matrix and the strictly
 *    lower triangular part of `a` is not referenced. With uplo='L' or 'l', the leading n by n lower
 *    triangular part of the array `a` must contain the lower triangular matrix and the strictly upper
 *    triangular part of `a` is not referenced. Note that when  diag='U' or 'u', the diagonal elements
 *    of a are not referenced either, but are assumed to be unity.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, n).
 *
 *    [out] x    - real-floating array of dimension (at least) (1+(n-1)*abs(incx)), the incremented array
 *    x must contain the n element right-hand side vector b, x is overwritten with the solution vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *              | 0 0 0 0 |
 *     a[4x4] = | 1 0 0 0 |
 *              | 2 3 0 0 |
 *              | 3 4 3 0 |
 *
 *     const real-floating a[] = {
 *          0, 0, 0, 0
 *        , 1, 0, 0, 0
 *        , 2, 3, 0, 0
 *        , 3, 4, 3, 0
 *     };
 *     real-floating x[] = { 1, 3, 11, 24 };
 *     mc_blas_?trsv('L', 'N', 'U', 4, a, 4, x, 1);
 *     on output -> x = { 1, 2, 3, 4 }
 *
 *              | 1 2 3 2 |
 *     a[4x4] = | 0 2 2 5 |
 *              | 0 0 3 3 |
 *              | 0 0 0 1 |
 *
 *     const real-floating a[] = {
 *          1, 2, 3, 2
 *        , 0, 2, 2, 5
 *        , 0, 0, 3, 3
 *        , 0, 0, 0, 1
 *     };
 *     real-floating x[] = { 5, 18, 32, 41 };
 *     mc_blas_?trsv('U', 'T', 'N', 4 , a, 4 , x, 1);
 *     on output -> x = { 5, 4, 3, 2 }
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
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_TRSV_H
#define MC_BLAS_TRSV_H

#pragma mark - mc_blas_strsv -

MC_TARGET_FUNC void mc_blas_strsv(const char uplo, const char trans, const char diag, const int n, const float * a, const int lda, float * x, const int incx)
{
	const float zero = 0.0f;

	float temp;
	int i, info, ix, j, jx, kx;
	int nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, n)) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	}
	if (info != 0) {
		mc_blas_xerbla("STRSV ", info);
		return;
	}

	if (n == 0) {
		return;
	}
	nounit = mc_blas_lsame(diag, 'N');

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, j) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = j - 1; i >= 1; --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
				}
			} else {
				jx = kx + (n - 1) * incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, jx) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, jx);
						ix   = jx;
						for (i = j - 1; i >= 1; --i) {
							ix = ix - incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
					jx = jx - incx;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = j + 1; i <= n; ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
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
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, jx);
						ix = jx;
						for (i = j + 1; i <= n; ++i) {
							ix                       = ix + incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
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
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, j) = temp;
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
					ix = kx;
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
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
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, j);
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_dtrsv -

MC_TARGET_FUNC void mc_blas_dtrsv(const char uplo, const char trans, const char diag, const int n, const double * a, const int lda, double * x, const int incx)
{
	const double zero = 0.0;

	double temp;
	int i, info, ix, j, jx, kx;
	int nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, n)) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	}
	if (info != 0) {
		mc_blas_xerbla("DTRSV ", info);
		return;
	}

	if (n == 0) {
		return;
	}
	nounit = mc_blas_lsame(diag, 'N');

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, j) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = j - 1; i >= 1; --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
				}
			} else {
				jx = kx + (n - 1) * incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, jx) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, jx);
						ix   = jx;
						for (i = j - 1; i >= 1; --i) {
							ix = ix - incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
					jx = jx - incx;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = j + 1; i <= n; ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
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
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, jx);
						ix = jx;
						for (i = j + 1; i <= n; ++i) {
							ix                       = ix + incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
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
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, j) = temp;
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
					ix = kx;
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
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
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, j);
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_ltrsv -

MC_TARGET_FUNC void mc_blas_ltrsv(const char uplo, const char trans, const char diag, const int n, const long double * a, const int lda, long double * x, const int incx)
{
	const long double zero = 0.0L;

	long double temp;
	int i, info, ix, j, jx, kx;
	int nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (lda < mc_maxmag(1, n)) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	}
	if (info != 0) {
		mc_blas_xerbla("LTRSV ", info);
		return;
	}

	if (n == 0) {
		return;
	}
	nounit = mc_blas_lsame(diag, 'N');

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, j) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = j - 1; i >= 1; --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
				}
			} else {
				jx = kx + (n - 1) * incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, jx) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, jx);
						ix   = jx;
						for (i = j - 1; i >= 1; --i) {
							ix = ix - incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
					jx = jx - incx;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = j + 1; i <= n; ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
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
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, j, j);
						}
						temp = mc_blas_vector_at(x, jx);
						ix = jx;
						for (i = j + 1; i <= n; ++i) {
							ix                       = ix + incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
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
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, j) = temp;
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
					ix = kx;
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
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
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, j);
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, j, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
				}
			}
		}
	}
}

#endif /* !MC_BLAS_TRSV_H */

/* EOF */