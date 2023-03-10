//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_syr2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?syr2 performs the symmetric rank 2 operation:
 *    a=alpha*x*y' + alpha*y*x' + a.
 *
 * \synopsis
 *    void ?syr2(uplo, n, alpha, x, incx, y, incy, a, lda)
 *    real-floating alpha
 *    int           incx, incy, lda, n
 *    char          uplo
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?syr2 performs the symmetric rank 2 operation: a=alpha*x*y' + alpha*y*x' + a where alpha
 *    is a scalar, `x` and `y` are n element vectors and `a` is an n by n symmetric matrix.
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
 *    [in]  y    - real-floating array of dimension (at least) (1+(n-1)*abs(incy)). The incremented array `y`
 *    must contain the vector `y`.
 *
 *    [in]  incy - int. Specifies the increment for the elements of `y`, incx must not be zero.
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

#ifndef MC_BLAS_SYR2_H
#define MC_BLAS_SYR2_H

#pragma mark - mc_blas_ssyr2 -

MC_TARGET_FUNC void mc_blas_ssyr2(const char uplo, const int n, const float alpha, const float * x, const int incx, const float * y, const int incy, float * a, const int lda)
{
	const float zero = 0.0f;

	float temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kx, ky;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, n)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("SSYR2 ", info);
		return;
	}

	if (n == 0 || alpha == zero) {
		return;
	}

	if (incx != 1 || incy != 1) {
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
		jx = kx;
		jy = ky;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero || mc_blas_vector_at(y, j) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, j);
					temp2 = alpha * mc_blas_vector_at(x, j);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
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
				if (mc_blas_vector_at(x, jx) != zero || mc_blas_vector_at(y, jy) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, jy);
					temp2 = alpha * mc_blas_vector_at(x, jx);
					ix    = kx;
					iy    = ky;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
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
				if (mc_blas_vector_at(x, j) != zero || mc_blas_vector_at(y, j) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, j);
					temp2 = alpha * mc_blas_vector_at(x, j);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
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
				if (mc_blas_vector_at(x, jx) != zero || mc_blas_vector_at(y, jy) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, jy);
					temp2 = alpha * mc_blas_vector_at(x, jx);
					ix    = jx;
					iy    = jy;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_dsyr2 -

MC_TARGET_FUNC void mc_blas_dsyr2(const char uplo, const int n, const double alpha, const double * x, const int incx, const double * y, const int incy, double * a, const int lda)
{
	const double zero = 0.0;

	double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kx, ky;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, n)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("DSYR2 ", info);
		return;
	}

	if (n == 0 || alpha == zero) {
		return;
	}

	if (incx != 1 || incy != 1) {
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
		jx = kx;
		jy = ky;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero || mc_blas_vector_at(y, j) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, j);
					temp2 = alpha * mc_blas_vector_at(x, j);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
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
				if (mc_blas_vector_at(x, jx) != zero || mc_blas_vector_at(y, jy) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, jy);
					temp2 = alpha * mc_blas_vector_at(x, jx);
					ix    = kx;
					iy    = ky;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
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
				if (mc_blas_vector_at(x, j) != zero || mc_blas_vector_at(y, j) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, j);
					temp2 = alpha * mc_blas_vector_at(x, j);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
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
				if (mc_blas_vector_at(x, jx) != zero || mc_blas_vector_at(y, jy) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, jy);
					temp2 = alpha * mc_blas_vector_at(x, jx);
					ix    = jx;
					iy    = jy;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_lsyr2 -

MC_TARGET_FUNC void mc_blas_lsyr2(const char uplo, const int n, const long double alpha, const long double * x, const int incx, const long double * y, const int incy, long double * a, const int lda)
{
	const long double zero = 0.0L;

	long double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kx, ky;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, n)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("LSYR2 ", info);
		return;
	}

	if (n == 0 || alpha == zero) {
		return;
	}

	if (incx != 1 || incy != 1) {
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
		jx = kx;
		jy = ky;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				if (mc_blas_vector_at(x, j) != zero || mc_blas_vector_at(y, j) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, j);
					temp2 = alpha * mc_blas_vector_at(x, j);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
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
				if (mc_blas_vector_at(x, jx) != zero || mc_blas_vector_at(y, jy) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, jy);
					temp2 = alpha * mc_blas_vector_at(x, jx);
					ix    = kx;
					iy    = ky;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
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
				if (mc_blas_vector_at(x, j) != zero || mc_blas_vector_at(y, j) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, j);
					temp2 = alpha * mc_blas_vector_at(x, j);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
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
				if (mc_blas_vector_at(x, jx) != zero || mc_blas_vector_at(y, jy) != zero) {
					temp1 = alpha * mc_blas_vector_at(y, jy);
					temp2 = alpha * mc_blas_vector_at(x, jx);
					ix    = jx;
					iy    = jy;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
}

/* \name
 *    ?syr2 performs the symmetric rank 2 operation:
 *    a=alpha*x*y_ + alpha*y*x_ + a.
 *
 * \synopsis
 *    void ?syr2(uplo, n, alpha, x, incx, y, incy, a, lda)
 *    complex alpha
 *    int     incx, incy, lda, n
 *    char    uplo
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?syr2 performs the symmetric rank 2 operation: a=alpha*x*y' + alpha*y*x' + a where alpha
 *    is a scalar, `x` and `y` are n element vectors and `a` is an n by n symmetric matrix.
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
 *    [in]  y    - complex array of dimension (at least) (1+(n-1)*abs(incy)). The incremented array `y`
 *    must contain the vector `y`.
 *
 *    [in]  incy - int. Specifies the increment for the elements of `y`, incx must not be zero.
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

#pragma mark - mc_blas_csyr2 -

MC_TARGET_FUNC void mc_blas_csyr2(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y, const int incy, mc_complex_float_t * a, const int lda)
{
	const mc_complex_float_t zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kx, ky;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, n)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("CSYR2 ", info);
		return;
	}

	if (n == 0 || mc_ciseqf(alpha, zero)) {
		return;
	}

	if (incx != 1 || incy != 1) {
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
		jx = kx;
		jy = ky;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_blas_vector_at(y, j);
				temp2 = mc_blas_vector_at(x, j);
				if (!mc_ciseqf(temp2, zero) || !mc_ciseqf(temp1, zero)) {
					temp1 = mc_cmulf(alpha, temp1);
					temp2 = mc_cmulf(alpha, temp2);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_caddf(mc_cmulf(mc_blas_vector_at(x, i), temp1), mc_cmulf(mc_blas_vector_at(y, i), temp2)));
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
				temp1 = mc_blas_vector_at(y, jy);
				temp2 = mc_blas_vector_at(x, jx);
				if (!mc_ciseqf(temp2, zero) || !mc_ciseqf(temp1, zero)) {
					temp1 = mc_cmulf(alpha, temp1);
					temp2 = mc_cmulf(alpha, temp2);
					ix    = kx;
					iy    = ky;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_caddf(mc_cmulf(mc_blas_vector_at(x, ix), temp1), mc_cmulf(mc_blas_vector_at(y, iy), temp2)));
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
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
				temp1 = mc_blas_vector_at(y, j);
				temp2 = mc_blas_vector_at(x, j);
				if (!mc_ciseqf(temp2, zero) || !mc_ciseqf(temp1, zero)) {
					temp1 = mc_cmulf(alpha, temp1);
					temp2 = mc_cmulf(alpha, temp2);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_caddf(mc_cmulf(mc_blas_vector_at(x, i), temp1), mc_cmulf(mc_blas_vector_at(y, i), temp2)));
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
				temp1 = mc_blas_vector_at(y, jy);
				temp2 = mc_blas_vector_at(x, jx);
				if (!mc_ciseqf(temp2, zero) || !mc_ciseqf(temp1, zero)) {
					temp1 = mc_cmulf(alpha, temp1);
					temp2 = mc_cmulf(alpha, temp2);
					ix    = jx;
					iy    = jy;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddf(mc_blas_matrix_at(a, lda, n, i, j), mc_caddf(mc_cmulf(mc_blas_vector_at(x, ix), temp1), mc_cmulf(mc_blas_vector_at(y, iy), temp2)));
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_zsyr2 -

MC_TARGET_FUNC void mc_blas_zsyr2(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y, const int incy, mc_complex_double_t * a, const int lda)
{
	const mc_complex_double_t zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kx, ky;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, n)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("ZSYR2 ", info);
		return;
	}

	if (n == 0 || mc_ciseq(alpha, zero)) {
		return;
	}

	if (incx != 1 || incy != 1) {
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
		jx = kx;
		jy = ky;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_blas_vector_at(y, j);
				temp2 = mc_blas_vector_at(x, j);
				if (!mc_ciseq(temp2, zero) || !mc_ciseq(temp1, zero)) {
					temp1 = mc_cmul(alpha, temp1);
					temp2 = mc_cmul(alpha, temp2);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cadd(mc_cmul(mc_blas_vector_at(x, i), temp1), mc_cmul(mc_blas_vector_at(y, i), temp2)));
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
				temp1 = mc_blas_vector_at(y, jy);
				temp2 = mc_blas_vector_at(x, jx);
				if (!mc_ciseq(temp2, zero) || !mc_ciseq(temp1, zero)) {
					temp1 = mc_cmul(alpha, temp1);
					temp2 = mc_cmul(alpha, temp2);
					ix    = kx;
					iy    = ky;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cadd(mc_cmul(mc_blas_vector_at(x, ix), temp1), mc_cmul(mc_blas_vector_at(y, iy), temp2)));
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
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
				temp1 = mc_blas_vector_at(y, j);
				temp2 = mc_blas_vector_at(x, j);
				if (!mc_ciseq(temp2, zero) || !mc_ciseq(temp1, zero)) {
					temp1 = mc_cmul(alpha, temp1);
					temp2 = mc_cmul(alpha, temp2);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cadd(mc_cmul(mc_blas_vector_at(x, i), temp1), mc_cmul(mc_blas_vector_at(y, i), temp2)));
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
				temp1 = mc_blas_vector_at(y, jy);
				temp2 = mc_blas_vector_at(x, jx);
				if (!mc_ciseq(temp2, zero) || !mc_ciseq(temp1, zero)) {
					temp1 = mc_cmul(alpha, temp1);
					temp2 = mc_cmul(alpha, temp2);
					ix    = jx;
					iy    = jy;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_cadd(mc_blas_matrix_at(a, lda, n, i, j), mc_cadd(mc_cmul(mc_blas_vector_at(x, ix), temp1), mc_cmul(mc_blas_vector_at(y, iy), temp2)));
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
}

#pragma mark - mc_blas_qsyr2 -

MC_TARGET_FUNC void mc_blas_qsyr2(const char uplo, const int n, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * x, const int incx, const mc_complex_long_double_t * y, const int incy, mc_complex_long_double_t * a, const int lda)
{
	const mc_complex_long_double_t zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, kx, ky;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	} else if (lda < mc_maxmag(1, n)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("QSYR2 ", info);
		return;
	}

	if (n == 0 || mc_ciseql(alpha, zero)) {
		return;
	}

	if (incx != 1 || incy != 1) {
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
		jx = kx;
		jy = ky;
	}

	if (mc_blas_lsame(uplo, 'U')) {
		if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (j = 1; j <= n; ++j) {
				temp1 = mc_blas_vector_at(y, j);
				temp2 = mc_blas_vector_at(x, j);
				if (!mc_ciseql(temp2, zero) || !mc_ciseql(temp1, zero)) {
					temp1 = mc_cmull(alpha, temp1);
					temp2 = mc_cmull(alpha, temp2);
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_caddl(mc_cmull(mc_blas_vector_at(x, i), temp1), mc_cmull(mc_blas_vector_at(y, i), temp2)));
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
				temp1 = mc_blas_vector_at(y, jy);
				temp2 = mc_blas_vector_at(x, jx);
				if (!mc_ciseql(temp2, zero) || !mc_ciseql(temp1, zero)) {
					temp1 = mc_cmull(alpha, temp1);
					temp2 = mc_cmull(alpha, temp2);
					ix    = kx;
					iy    = ky;
					for (i = 1; i <= j; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_caddl(mc_cmull(mc_blas_vector_at(x, ix), temp1), mc_cmull(mc_blas_vector_at(y, iy), temp2)));
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
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
				temp1 = mc_blas_vector_at(y, j);
				temp2 = mc_blas_vector_at(x, j);
				if (!mc_ciseql(temp2, zero) || !mc_ciseql(temp1, zero)) {
					temp1 = mc_cmull(alpha, temp1);
					temp2 = mc_cmull(alpha, temp2);
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_caddl(mc_cmull(mc_blas_vector_at(x, i), temp1), mc_cmull(mc_blas_vector_at(y, i), temp2)));
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
				temp1 = mc_blas_vector_at(y, jy);
				temp2 = mc_blas_vector_at(x, jx);
				if (!mc_ciseql(temp2, zero) || !mc_ciseql(temp1, zero)) {
					temp1 = mc_cmull(alpha, temp1);
					temp2 = mc_cmull(alpha, temp2);
					ix    = jx;
					iy    = jy;
					for (i = j; i <= n; ++i) {
						mc_blas_matrix_at(a, lda, n, i, j) = mc_caddl(mc_blas_matrix_at(a, lda, n, i, j), mc_caddl(mc_cmull(mc_blas_vector_at(x, ix), temp1), mc_cmull(mc_blas_vector_at(y, iy), temp2)));
						ix                                 = ix + incx;
						iy                                 = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
}

#endif /* !MC_BLAS_SYR2_H */

/* EOF */