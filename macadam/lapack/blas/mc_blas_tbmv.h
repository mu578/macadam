//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_tbmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?tbmv performs the matrix-vector operation:
 *    x=a*x or x=a'*x.
 *
 * \synopsis
 *    void ?tbmv(uplo, trans, diag, n, k, a, lda, x, incx)
 *    int           incx, k, lda, n
 *    char          diag, trans, uplo
 *    real-floating a(lda, *), x(*)
 *
 * \purpose
 *    ?tbmv performs the matrix-vector operation: x=a*x or x=a'*x where `x` is an n element vector and
 *    `a` is an n by n unit or non-unit, upper or lower triangular band matrix, with (k + 1) diagonals.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the matrix is an upper or lower triangular matrix as follows:
 *    uplo='U' or 'u', `a` is an upper triangular matrix.
 *    uplo='L' or 'l', `a` is a lower triangular matrix.
 *
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' x=a*x.
 *    trans='T' or 't' x=a'*x.
 *    trans='C' or 'c' x=a'*x.
 *
 *    [in] diag  - char. Specifies whether or not `a` is unit triangular as follows:
 *    diag='U' or 'u' `a` is assumed to be unit triangular.
 *    diag='N' or 'n' `a` is not assumed to be unit triangular.
 *
 *    [in] n     - int. Specifies the order of the matrix `a`, n must be at least zero.
 *    [in] k     - int. With uplo='U' or 'u', k specifies the number of super-diagonals of the matrix `a`.
 *    i.e the half-bandwidth of matrix `a`. With uplo='L' or 'l', k specifies the number of sub-diagonals
 *    of the matrix `a`.k must satisfy 0 < k.
 *
 *    [in] a     - real-floating array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column-by-column with the main diagonal of the matrix in row (k) of the
 *    array, the first super-diagonal starting at position 1 in row (k - 1), and so on. Elements in the array `a` which do
 *    not correspond to elements in the triangular band matrix (such as the top left k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the upper triangular band part of the matrix of
 *    coefficients. The matrix must be supplied row-by-row, with the leading diagonal of the matrix in column 0 of the array,
 *    the first super-diagonal starting at position 0 in column 1, and so on. The bottom right k by k triangle of array `a`
 *    is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column by column with the main diagonal of the matrix in row 0 of the
 *    array, the first sub-diagonal starting at position 0 in row 1, and so on. Elements in the array `a` that do not
 *    correspond to elements in the triangular band matrix (such as the bottom right k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the lower triangular band part of the matrix of
 *    coefficients, supplied row by row, with the leading diagonal of the matrix in column k of the array, the first sub-
 *    diagonal starting at position 1 in column (k - 1), and so on. The top left k by k triangle of array `a` is not
 *    referenced.
 *
 *    @see: \examples section about Triangular Band Matrix storage.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least (k + 1).
 *
 *    [out] x    - real-floating array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must contain the
 *    n element vector `x`. `x` is overwritten with the transformed vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *     - A triangular band matrix `a` of n rows and n columns with k sub/super-diagonals and leading dimension lda. @note: using
 *       @c-layout, be aware that the true-memory leading or split dimension becomes the number of columns (The leading row syntax
 *       is only kept for the understanding of the transposed applied operations in reference to Fortran). For `a` a m by n matrix,
 *       the rational subscript operation should be as the following:
 *
 *       int row, col;
 *       for (row = 0; row < m; row++) {
 *          for (col = 0; col < n; col++) {
 *             a[(n * row) + col] = ...;
 *          }
 *       }
 *
 *     - Given uplo='U' or 'u': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       real-floating b[ldb * n] = { ... };
 *       real-floating a[lda * n] = { 0 };
 *
 *       int i, j, m, k;
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = k + 1 - j
 *          for (i = mc_maxmag(1, j - k); i <= j; ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = k - j;
 *          for (i = mc_maxmag(0, j - k); i <= j; i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = -i;
 *          for (j = i; j < mc_minmag(n, i + k + 1); j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 *     - Given uplo='L' or 'l': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       real-floating b[ldb * n] = { ... };
 *       real-floating a[lda * n] = { 0 };
 *
 *       int i, j, m, k;
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = 1 - j
 *          for (i = j; i <= mc_minmag(n, j + k); ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = -j;
 *          for (i = j; i < mc_minmag(n, j + k + 1); i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = k - i;
 *          for (j = mc_maxmag(0, i - k); j <= i; j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
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

#ifndef MC_BLAS_TBMV_H
#define MC_BLAS_TBMV_H

#pragma mark - mc_blas_stbmv -

MC_TARGET_FUNC void mc_blas_stbmv(const char uplo, const char trans, const char diag, const int n, const int k, const float * a, const int lda, float * x, const int incx)
{
	const float zero = 0.0f;

	float temp;
	int i, info, ix, j, jx, kplus1, kx, l;
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
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("STBMV ", info);
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
			kplus1 = k + 1;
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
						temp = mc_blas_vector_at(x, j);
						l    = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, kplus1, j);
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
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						l    = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
					}
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
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
					if (mc_blas_vector_at(x, j) != zero) {
						temp = mc_blas_vector_at(x, j);
						l    = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, 1, j);
						}
					}
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						l    = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, 1, j);
						}
					}
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
					l    = kplus1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
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
					kx   = kx - incx;
					ix   = kx;
					l    = kplus1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					l    = 1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, 1, j);
					}
					for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
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
					kx   = kx + incx;
					ix   = kx;
					l    = 1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, 1, j);
					}
					for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_dtbmv -

MC_TARGET_FUNC void mc_blas_dtbmv(const char uplo, const char trans, const char diag, const int n, const int k, const double * a, const int lda, double * x, const int incx)
{
	const double zero = 0.0;

	double temp;
	int i, info, ix, j, jx, kplus1, kx, l;
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
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("DTBMV ", info);
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
			kplus1 = k + 1;
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
						temp = mc_blas_vector_at(x, j);
						l    = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, kplus1, j);
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
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						l    = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
					}
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
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
					if (mc_blas_vector_at(x, j) != zero) {
						temp = mc_blas_vector_at(x, j);
						l    = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, 1, j);
						}
					}
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						l    = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, 1, j);
						}
					}
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
					l    = kplus1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
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
					kx   = kx - incx;
					ix   = kx;
					l    = kplus1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					l    = 1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, 1, j);
					}
					for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
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
					kx   = kx + incx;
					ix   = kx;
					l    = 1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, 1, j);
					}
					for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_ltbmv -

MC_TARGET_FUNC void mc_blas_ltbmv(const char uplo, const char trans, const char diag, const int n, const int k, const long double * a, const int lda, long double * x, const int incx)
{
	const long double zero = 0.0;

	long double temp;
	int i, info, ix, j, jx, kplus1, kx, l;
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
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("LTBMV ", info);
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
			kplus1 = k + 1;
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
						temp = mc_blas_vector_at(x, j);
						l    = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, kplus1, j);
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
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						l    = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
					}
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
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
					if (mc_blas_vector_at(x, j) != zero) {
						temp = mc_blas_vector_at(x, j);
						l    = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, 1, j);
						}
					}
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						l    = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, 1, j);
						}
					}
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
					l    = kplus1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
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
					kx   = kx - incx;
					ix   = kx;
					l    = kplus1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					l    = 1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, 1, j);
					}
					for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
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
					kx   = kx + incx;
					ix   = kx;
					l    = 1 - j;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, 1, j);
					}
					for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

/* \name
 *    ?tbmv performs the matrix-vector operation:
 *    x=a*x or x=a'*x or x=a_*x.
 *
 * \synopsis
 *    void ?tbmv(uplo, trans, diag, n, k, a, lda, x, incx)
 *    int     incx, k, lda, n
 *    char    diag, trans, uplo
 *    complex a(lda, *), x(*)
 *
 * \purpose
 *    ?tbmv performs the matrix-vector operation: x=a*x or x=a'*x or x=a_*x where `x` is an n element vector
 *    and `a` is an n by n unit or non-unit, upper or lower triangular band matrix, with (k + 1) diagonals.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the matrix is an upper or lower triangular matrix as follows:
 *    uplo='U' or 'u', `a` is an upper triangular matrix.
 *    uplo='L' or 'l', `a` is a lower triangular matrix.
 *
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' x=a*x.
 *    trans='T' or 't' x=a'*x.
 *    trans='C' or 'c' x=a_*x.
 *
 *    [in] diag  - char. Specifies whether or not `a` is unit triangular as follows:
 *    diag='U' or 'u' `a` is assumed to be unit triangular.
 *    diag='N' or 'n' `a` is not assumed to be unit triangular.
 *
 *    [in] n     - int. Specifies the order of the matrix `a`, n must be at least zero.
 *    [in] k     - int. With uplo='U' or 'u', k specifies the number of super-diagonals of the matrix `a`.
 *    i.e the half-bandwidth of matrix `a`. With uplo='L' or 'l', k specifies the number of sub-diagonals
 *    of the matrix `a`.k must satisfy 0 < k.
 *
 *    [in] a     - complex array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column-by-column with the main diagonal of the matrix in row (k) of the
 *    array, the first super-diagonal starting at position 1 in row (k - 1), and so on. Elements in the array `a` which do
 *    not correspond to elements in the triangular band matrix (such as the top left k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the upper triangular band part of the matrix of
 *    coefficients. The matrix must be supplied row-by-row, with the leading diagonal of the matrix in column 0 of the array,
 *    the first super-diagonal starting at position 0 in column 1, and so on. The bottom right k by k triangle of array `a`
 *    is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column by column with the main diagonal of the matrix in row 0 of the
 *    array, the first sub-diagonal starting at position 0 in row 1, and so on. Elements in the array `a` that do not
 *    correspond to elements in the triangular band matrix (such as the bottom right k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the lower triangular band part of the matrix of
 *    coefficients, supplied row by row, with the leading diagonal of the matrix in column k of the array, the first sub-
 *    diagonal starting at position 1 in column (k - 1), and so on. The top left k by k triangle of array `a` is not
 *    referenced.
 *
 *    @see: \examples section about Triangular Band Matrix storage.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least (k + 1).
 *
 *    [out] x    - complex array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must contain the
 *    n element vector `x`. `x` is overwritten with the transformed vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *     - A triangular band matrix `a` of n rows and n columns with k sub/super-diagonals and leading dimension lda. @note: using
 *       @c-layout, be aware that the true-memory leading or split dimension becomes the number of columns (The leading row syntax
 *       is only kept for the understanding of the transposed applied operations in reference to Fortran). For `a` a m by n matrix,
 *       the rational subscript operation should be as the following:
 *
 *       int row, col;
 *       for (row = 0; row < m; row++) {
 *          for (col = 0; col < n; col++) {
 *             a[(n * row) + col] = ...;
 *          }
 *       }
 *
 *     - Given uplo='U' or 'u': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       complex b[ldb * n] = { ... };
 *       complex a[lda * n] = { 0 };
 *
 *       int i, j, m, k = lda - 1; // must statisfy ldb >= lda >= n > k.
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = k + 1 - j
 *          for (i = mc_maxmag(1, j - k); i <= j; ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = k - j;
 *          for (i = mc_maxmag(0, j - k); i <= j; i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = -i;
 *          for (j = i; j < mc_minmag(n, i + k + 1); j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 *     - Given uplo='L' or 'l': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       complex b[ldb * n] = { ... };
 *       complex a[lda * n] = { 0 };
 *
 *       int i, j, m, k = lda - 1; // must statisfy ldb >= lda >= n > k.
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = 1 - j
 *          for (i = j; i <= mc_minmag(n, j + k); ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = -j;
 *          for (i = j; i < mc_minmag(n, j + k + 1); i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = k - i;
 *          for (j = mc_maxmag(0, i - k); j <= i; j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
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

#pragma mark - mc_blas_ctbmv -

MC_TARGET_FUNC void mc_blas_ctbmv(const char uplo, const char trans, const char diag, const int n, const int k, const mc_complex_float_t * a, const int lda, mc_complex_float_t * x, const int incx)
{
	const mc_complex_float_t zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp;
	int i, info, ix, j, jx, kplus1, kx, l;
	int noconj, nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("CTBMV ", info);
		return;
	}

	if (n == 0) {
		return;
	}

	noconj = mc_blas_lsame(trans, 'T');
	nounit = mc_blas_lsame(diag, 'N');
	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
						l = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_caddf(mc_blas_vector_at(x, i), mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_cmulf(mc_blas_vector_at(x, j), mc_blas_matrix_at(a, lda, n, kplus1, j));
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
						ix = kx;
						l  = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_caddf(mc_blas_vector_at(x, ix), mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_cmulf(mc_blas_vector_at(x, jx), mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
					}
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
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
					if (!mc_ciseqf(temp, zero)) {
						l = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_caddf(mc_blas_vector_at(x, i), mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_cmulf(mc_blas_vector_at(x, j), mc_blas_matrix_at(a, lda, n, 1, j));
						}
					}
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
					if (!mc_ciseqf(temp, zero)) {
						ix = kx;
						l  = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_caddf(mc_blas_vector_at(x, ix), mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_cmulf(mc_blas_vector_at(x, jx), mc_blas_matrix_at(a, lda, n, 1, j));
						}
					}
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
					l    = kplus1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
						}
					} else {
						if (nounit) {
							temp = mc_cmulf(temp, mc_conjf(mc_blas_matrix_at(a, lda, n, kplus1, j)));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, i)));
						}
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
					kx   = kx - incx;
					ix   = kx;
					l    = kplus1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
							ix   = ix - incx;
						}
					} else {
						if (nounit) {
							temp = mc_cmulf(temp, mc_conjf(mc_blas_matrix_at(a, lda, n, kplus1, j)));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, ix)));
							ix   = ix - incx;
						}
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					l    = 1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, 1, j));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
						}
					} else {
						if (nounit) {
							temp = mc_cmulf(temp, mc_conjf(mc_blas_matrix_at(a, lda, n, 1, j)));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, i)));
						}
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
					kx   = kx + incx;
					ix   = kx;
					l    = 1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmulf(temp, mc_blas_matrix_at(a, lda, n, 1, j));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddf(temp, mc_cmulf(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
							ix   = ix + incx;
						}
					} else {
						if (nounit) {
							temp = mc_cmulf(temp, mc_conjf(mc_blas_matrix_at(a, lda, n, 1, j)));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, ix)));
							ix   = ix + incx;
						}
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_ztbmv -

MC_TARGET_FUNC void mc_blas_ztbmv(const char uplo, const char trans, const char diag, const int n, const int k, const mc_complex_double_t * a, const int lda, mc_complex_double_t * x, const int incx)
{
	const mc_complex_double_t zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp;
	int i, info, ix, j, jx, kplus1, kx, l;
	int noconj, nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("ZTBMV ", info);
		return;
	}

	if (n == 0) {
		return;
	}

	noconj = mc_blas_lsame(trans, 'T');
	nounit = mc_blas_lsame(diag, 'N');
	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
						l = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_cadd(mc_blas_vector_at(x, i), mc_cmul(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_cmul(mc_blas_vector_at(x, j), mc_blas_matrix_at(a, lda, n, kplus1, j));
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
						ix = kx;
						l  = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_cadd(mc_blas_vector_at(x, ix), mc_cmul(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_cmul(mc_blas_vector_at(x, jx), mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
					}
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
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
					if (!mc_ciseq(temp, zero)) {
						l = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_cadd(mc_blas_vector_at(x, i), mc_cmul(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_cmul(mc_blas_vector_at(x, j), mc_blas_matrix_at(a, lda, n, 1, j));
						}
					}
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
					if (!mc_ciseq(temp, zero)) {
						ix = kx;
						l  = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_cadd(mc_blas_vector_at(x, ix), mc_cmul(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_cmul(mc_blas_vector_at(x, jx), mc_blas_matrix_at(a, lda, n, 1, j));
						}
					}
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
					l    = kplus1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmul(temp, mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
						}
					} else {
						if (nounit) {
							temp = mc_cmul(temp, mc_conj(mc_blas_matrix_at(a, lda, n, kplus1, j)));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, i)));
						}
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
					kx   = kx - incx;
					ix   = kx;
					l    = kplus1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmul(temp, mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
							ix   = ix - incx;
						}
					} else {
						if (nounit) {
							temp = mc_cmul(temp, mc_conj(mc_blas_matrix_at(a, lda, n, kplus1, j)));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, ix)));
							ix   = ix - incx;
						}
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					l    = 1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmul(temp, mc_blas_matrix_at(a, lda, n, 1, j));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
						}
					} else {
						if (nounit) {
							temp = mc_cmul(temp, mc_conj(mc_blas_matrix_at(a, lda, n, 1, j)));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, i)));
						}
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
					kx   = kx + incx;
					ix   = kx;
					l    = 1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmul(temp, mc_blas_matrix_at(a, lda, n, 1, j));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_cadd(temp, mc_cmul(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
							ix   = ix + incx;
						}
					} else {
						if (nounit) {
							temp = mc_cmul(temp, mc_conj(mc_blas_matrix_at(a, lda, n, 1, j)));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, ix)));
							ix   = ix + incx;
						}
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_qtbmv -

MC_TARGET_FUNC void mc_blas_qtbmv(const char uplo, const char trans, const char diag, const int n, const int k, const mc_complex_long_double_t * a, const int lda, mc_complex_long_double_t * x, const int incx)
{
	const mc_complex_long_double_t zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp;
	int i, info, ix, j, jx, kplus1, kx, l;
	int noconj, nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("QTBMV ", info);
		return;
	}

	if (n == 0) {
		return;
	}

	noconj = mc_blas_lsame(trans, 'T');
	nounit = mc_blas_lsame(diag, 'N');
	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
						l = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_caddl(mc_blas_vector_at(x, i), mc_cmull(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_cmull(mc_blas_vector_at(x, j), mc_blas_matrix_at(a, lda, n, kplus1, j));
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
						ix = kx;
						l  = kplus1 - j;
						for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_caddl(mc_blas_vector_at(x, ix), mc_cmull(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_cmull(mc_blas_vector_at(x, jx), mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
					}
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
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
					if (!mc_ciseql(temp, zero)) {
						l = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_caddl(mc_blas_vector_at(x, i), mc_cmull(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_cmull(mc_blas_vector_at(x, j), mc_blas_matrix_at(a, lda, n, 1, j));
						}
					}
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
					if (!mc_ciseql(temp, zero)) {
						ix = kx;
						l  = 1 - j;
						for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_caddl(mc_blas_vector_at(x, ix), mc_cmull(temp, mc_blas_matrix_at(a, lda, n, l + i, j)));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_cmull(mc_blas_vector_at(x, jx), mc_blas_matrix_at(a, lda, n, 1, j));
						}
					}
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
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
					l    = kplus1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmull(temp, mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
						}
					} else {
						if (nounit) {
							temp = mc_cmull(temp, mc_conjl(mc_blas_matrix_at(a, lda, n, kplus1, j)));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, i)));
						}
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
					kx   = kx - incx;
					ix   = kx;
					l    = kplus1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmull(temp, mc_blas_matrix_at(a, lda, n, kplus1, j));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
							ix   = ix - incx;
						}
					} else {
						if (nounit) {
							temp = mc_cmull(temp, mc_conjl(mc_blas_matrix_at(a, lda, n, kplus1, j)));
						}
						for (i = j - 1; i >= mc_maxmag(1, j - k); --i) {
							temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, ix)));
							ix   = ix - incx;
						}
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					l    = 1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmull(temp, mc_blas_matrix_at(a, lda, n, 1, j));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, i)));
						}
					} else {
						if (nounit) {
							temp = mc_cmull(temp, mc_conjl(mc_blas_matrix_at(a, lda, n, 1, j)));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, i)));
						}
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
					kx   = kx + incx;
					ix   = kx;
					l    = 1 - j;
					if (noconj) {
						if (nounit) {
							temp = mc_cmull(temp, mc_blas_matrix_at(a, lda, n, 1, j));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddl(temp, mc_cmull(mc_blas_matrix_at(a, lda, n, l + i, j), mc_blas_vector_at(x, ix)));
							ix   = ix + incx;
						}
					} else {
						if (nounit) {
							temp = mc_cmull(temp, mc_conjl(mc_blas_matrix_at(a, lda, n, 1, j)));
						}
						for (i = j + 1; i <= mc_minmag(n, j + k); ++i) {
							temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_matrix_at(a, lda, n, l + i, j)), mc_blas_vector_at(x, ix)));
							ix   = ix + incx;
						}
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#endif /* !MC_BLAS_TBMV_H */

/* EOF */