//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_spmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?spmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?spmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
 *    real-floating alpha, beta
 *    int           incx, incy, n
 *    char          uplo
 *    real-floating ap(*), x(*), y(*)
 *
 * \purpose
 *    ?spmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `ap` is an n by n symmetric matrix, supplied in
 *    packed form of dimension (at least) ((n*(n+1))/2).
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the matrix `a` is
 *    supplied in the packed array `ap` as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` supplied in `ap`.
 *    uplo='L' or 'l', the lower triangular part of `a` supplied in `ap`.
 *
 *    [in] n     - int. Specifies the order of the symmetric matrix `a`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] ap    - real-floating array of dimension (at least) ((n*(n+1))/2).
 *    With uplo='U' or 'u', the array `ap` must contain the upper triangular part of the symmetric matrix
 *    packed sequentially, if column-major layout: column by column, so that ap(1) contains a(1,1), ap(2)
 *    and ap(3) contain a(1,2) and a(2,2) respectively, and so on, else row by row following the same logic.
 *
 *    With uplo='L' or 'l', the array `a`p must contain the lower triangular part of the symmetric matrix
 *    packed sequentially, if column-major layout: column by column, so that ap(1) contains a(1,1), ap(2)
 *    and ap(3) contain a(2,1)  and a(3,1) respectively, and so on, else row by row following the same logic.
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

#ifndef MC_BLAS_SPMV_H
#define MC_BLAS_SPMV_H

#pragma mark - mc_blas_sspmv -

MC_TARGET_FUNC void mc_blas_sspmv(const char uplo, const int n, const float alpha, const float * ap, const float * x, const int incx, const float beta, float * y, const int incy)
{
	const float one = 1.0f, zero = 0.0f;

	float temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_ap = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_ap = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_ap, 'U') && !mc_blas_lsame(uplo_ap, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 6;
	} else if (incy == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("SSPMV ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_ap, 'U')) {
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
				k     = kk;
				for (i = 1; i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                   = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + temp1 * mc_blas_vector_at(ap, kk + j - 1) + alpha * temp2;
				kk                      = kk + j;
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
				for (k = kk; k <= (kk + j - 2); ++k) {
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                    = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + temp1 * mc_blas_vector_at(ap, kk + j - 1) + alpha * temp2;
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + j;
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
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (temp1 * mc_blas_vector_at(ap, kk));
				k                       = kk + 1;
				for (i = j + 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                   = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (alpha * temp2);
				kk                      = kk + (n - j + 1);
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
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (temp1 * mc_blas_vector_at(ap, kk));
				ix                       = jx;
				iy                       = jy;
				for (k = (kk + 1); k <= (kk + n - j); ++k) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                    = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (alpha * temp2);
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + (n - j + 1);
			}
		}
	}
}

#pragma mark - mc_blas_dspmv -

MC_TARGET_FUNC void mc_blas_dspmv(const char uplo, const int n, const double alpha, const double * ap, const double * x, const int incx, const double beta, double * y, const int incy)
{
	const double one = 1.0, zero = 0.0;

	double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_ap = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_ap = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_ap, 'U') && !mc_blas_lsame(uplo_ap, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 6;
	} else if (incy == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("DSPMV ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_ap, 'U')) {
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
				k     = kk;
				for (i = 1; i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                   = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + temp1 * mc_blas_vector_at(ap, kk + j - 1) + alpha * temp2;
				kk                      = kk + j;
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
				for (k = kk; k <= (kk + j - 2); ++k) {
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                    = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + temp1 * mc_blas_vector_at(ap, kk + j - 1) + alpha * temp2;
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + j;
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
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (temp1 * mc_blas_vector_at(ap, kk));
				k                       = kk + 1;
				for (i = j + 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                   = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (alpha * temp2);
				kk                      = kk + (n - j + 1);
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
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (temp1 * mc_blas_vector_at(ap, kk));
				ix                       = jx;
				iy                       = jy;
				for (k = (kk + 1); k <= (kk + n - j); ++k) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                    = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (alpha * temp2);
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + (n - j + 1);
			}
		}
	}
}

#pragma mark - mc_blas_lspmv -

MC_TARGET_FUNC void mc_blas_lspmv(const char uplo, const int n, const long double alpha, const long double * ap, const long double * x, const int incx, const long double beta, long double * y, const int incy)
{
	const long double one = 1.0L, zero = 0.0L;

	long double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_ap = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_ap = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_ap, 'U') && !mc_blas_lsame(uplo_ap, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 6;
	} else if (incy == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("DSPMV ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_ap, 'U')) {
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
				k     = kk;
				for (i = 1; i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                   = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + temp1 * mc_blas_vector_at(ap, kk + j - 1) + alpha * temp2;
				kk                      = kk + j;
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
				for (k = kk; k <= (kk + j - 2); ++k) {
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                    = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + temp1 * mc_blas_vector_at(ap, kk + j - 1) + alpha * temp2;
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + j;
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
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (temp1 * mc_blas_vector_at(ap, kk));
				k                       = kk + 1;
				for (i = j + 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_blas_vector_at(y, i) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                   = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_blas_vector_at(y, j) + (alpha * temp2);
				kk                      = kk + (n - j + 1);
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
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (temp1 * mc_blas_vector_at(ap, kk));
				ix                       = jx;
				iy                       = jy;
				for (k = (kk + 1); k <= (kk + n - j); ++k) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_blas_vector_at(y, iy) + (temp1 * mc_blas_vector_at(ap, k));
					temp2                    = temp2 + (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
				}
				mc_blas_vector_at(y, jy) = mc_blas_vector_at(y, jy) + (alpha * temp2);
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + (n - j + 1);
			}
		}
	}
}

/* \name
 *    ?spmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?spmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
 *    complex alpha, beta
 *    int     incx, incy, n
 *    char    uplo
 *    complex ap(*), x(*), y(*)
 *
 * \purpose
 *    ?spmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `ap` is an n by n symmetric matrix, supplied in
 *    packed form of dimension (at least) ((n*(n+1))/2).
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the matrix `a` is
 *    supplied in the packed array `ap` as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` supplied in `ap`.
 *    uplo='L' or 'l', the lower triangular part of `a` supplied in `ap`.
 *
 *    [in] n     - int. Specifies the order of the symmetric matrix `a`, n must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] ap    - complex array of dimension (at least) ((n*(n+1))/2).
 *    With uplo='U' or 'u', the array `ap` must contain the upper triangular part of the symmetric matrix
 *    packed sequentially, if column-major layout: column by column, so that ap(1) contains a(1,1), ap(2)
 *    and ap(3) contain a(1,2) and a(2,2) respectively, and so on, else row by row following the same logic.
 *
 *    With uplo='L' or 'l', the array `a`p must contain the lower triangular part of the symmetric matrix
 *    packed sequentially, if column-major layout: column by column, so that ap(1) contains a(1,1), ap(2)
 *    and ap(3) contain a(2,1)  and a(3,1) respectively, and so on, else row by row following the same logic.
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

#pragma mark - mc_blas_cspmv -

MC_TARGET_FUNC void mc_blas_cspmv(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * ap, const mc_complex_float_t * x, const int incx, const mc_complex_float_t beta, mc_complex_float_t * y, const int incy)
{
	const mc_complex_float_t one = mc_cmplxf(1.0f, 0.0f), zero = mc_cmplxf(0.0f, 0.0f);

	mc_complex_float_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_ap = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_ap = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_ap, 'U') && !mc_blas_lsame(uplo_ap, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 6;
	} else if (incy == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("CSPMV ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_ap, 'U')) {
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
				k     = kk;
				for (i = 1; i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_caddf(mc_blas_vector_at(y, i), mc_cmulf(temp1, mc_blas_vector_at(ap, k)));
					temp2                   = mc_caddf(temp2, mc_cmulf(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, i)));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_caddf(mc_blas_vector_at(y, j), mc_caddf(mc_cmulf(temp1, mc_blas_vector_at(ap, kk + j - 1)), mc_cmulf(alpha, temp2)));
				kk                      = kk + j;
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
				for (k = kk; k <= (kk + j - 2); ++k) {
					mc_blas_vector_at(y, iy) = mc_caddf(mc_blas_vector_at(y, iy), mc_cmulf(temp1, mc_blas_vector_at(ap, k)));
					temp2                    = mc_caddf(temp2, mc_cmulf(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, ix)));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_caddf(mc_blas_vector_at(y, jy), mc_caddf(mc_cmulf(temp1, mc_blas_vector_at(ap, kk + j - 1)), mc_cmulf(alpha, temp2)));
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + j;
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
				mc_blas_vector_at(y, j) = mc_caddf(mc_blas_vector_at(y, j), mc_cmulf(temp1, mc_blas_vector_at(ap, kk)));
				k                       = kk + 1;
				for (i = j + 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_caddf(mc_blas_vector_at(y, i), mc_cmulf(temp1, mc_blas_vector_at(ap, k)));
					temp2                   = mc_caddf(temp2, mc_cmulf(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, i)));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_caddf(mc_blas_vector_at(y, j), mc_cmulf(alpha, temp2));
				kk                      = kk + (n - j + 1);
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
				mc_blas_vector_at(y, jy) = mc_caddf(mc_blas_vector_at(y, jy), mc_cmulf(temp1, mc_blas_vector_at(ap, kk)));
				ix                       = jx;
				iy                       = jy;
				for (k = (kk + 1); k <= (kk + n - j); ++k) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_caddf(mc_blas_vector_at(y, iy), mc_cmulf(temp1, mc_blas_vector_at(ap, k)));
					temp2                    = mc_caddf(temp2, mc_cmulf(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, ix)));
				}
				mc_blas_vector_at(y, jy) = mc_caddf(mc_blas_vector_at(y, jy), mc_cmulf(alpha, temp2));
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + (n - j + 1);
			}
		}
	}
}

#pragma mark - mc_blas_zspmv -

MC_TARGET_FUNC void mc_blas_zspmv(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * ap, const mc_complex_double_t * x, const int incx, const mc_complex_double_t beta, mc_complex_double_t * y, const int incy)
{
	const mc_complex_double_t one = mc_cmplx(1.0, 0.0), zero = mc_cmplx(0.0, 0.0);

	mc_complex_double_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_ap = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_ap = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_ap, 'U') && !mc_blas_lsame(uplo_ap, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 6;
	} else if (incy == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("ZSPMV ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_ap, 'U')) {
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
				k     = kk;
				for (i = 1; i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_cadd(mc_blas_vector_at(y, i), mc_cmul(temp1, mc_blas_vector_at(ap, k)));
					temp2                   = mc_cadd(temp2, mc_cmul(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, i)));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_cadd(mc_blas_vector_at(y, j), mc_cadd(mc_cmul(temp1, mc_blas_vector_at(ap, kk + j - 1)), mc_cmul(alpha, temp2)));
				kk                      = kk + j;
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
				for (k = kk; k <= (kk + j - 2); ++k) {
					mc_blas_vector_at(y, iy) = mc_cadd(mc_blas_vector_at(y, iy), mc_cmul(temp1, mc_blas_vector_at(ap, k)));
					temp2                    = mc_cadd(temp2, mc_cmul(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, ix)));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_cadd(mc_blas_vector_at(y, jy), mc_cadd(mc_cmul(temp1, mc_blas_vector_at(ap, kk + j - 1)), mc_cmul(alpha, temp2)));
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + j;
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
				mc_blas_vector_at(y, j) = mc_cadd(mc_blas_vector_at(y, j), mc_cmul(temp1, mc_blas_vector_at(ap, kk)));
				k                       = kk + 1;
				for (i = j + 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_cadd(mc_blas_vector_at(y, i), mc_cmul(temp1, mc_blas_vector_at(ap, k)));
					temp2                   = mc_cadd(temp2, mc_cmul(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, i)));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_cadd(mc_blas_vector_at(y, j), mc_cmul(alpha, temp2));
				kk                      = kk + (n - j + 1);
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
				mc_blas_vector_at(y, jy) = mc_cadd(mc_blas_vector_at(y, jy), mc_cmul(temp1, mc_blas_vector_at(ap, kk)));
				ix                       = jx;
				iy                       = jy;
				for (k = (kk + 1); k <= (kk + n - j); ++k) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_cadd(mc_blas_vector_at(y, iy), mc_cmul(temp1, mc_blas_vector_at(ap, k)));
					temp2                    = mc_cadd(temp2, mc_cmul(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, ix)));
				}
				mc_blas_vector_at(y, jy) = mc_cadd(mc_blas_vector_at(y, jy), mc_cmul(alpha, temp2));
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + (n - j + 1);
			}
		}
	}
}

#pragma mark - mc_blas_qspmv -

MC_TARGET_FUNC void mc_blas_qspmv(const char uplo, const int n, const mc_complex_long_double_t alpha, const mc_complex_long_double_t * ap, const mc_complex_long_double_t * x, const int incx, const mc_complex_long_double_t beta, mc_complex_long_double_t * y, const int incy)
{
	const mc_complex_long_double_t one = mc_cmplxl(1.0L, 0.0L), zero = mc_cmplxl(0.0L, 0.0L);

	mc_complex_long_double_t temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_ap = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_ap = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_ap, 'U') && !mc_blas_lsame(uplo_ap, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 6;
	} else if (incy == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("QSPMV ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_ap, 'U')) {
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
				k     = kk;
				for (i = 1; i <= (j - 1); ++i) {
					mc_blas_vector_at(y, i) = mc_caddl(mc_blas_vector_at(y, i), mc_cmull(temp1, mc_blas_vector_at(ap, k)));
					temp2                   = mc_caddl(temp2, mc_cmull(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, i)));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_caddl(mc_blas_vector_at(y, j), mc_caddl(mc_cmull(temp1, mc_blas_vector_at(ap, kk + j - 1)), mc_cmull(alpha, temp2)));
				kk                      = kk + j;
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
				for (k = kk; k <= (kk + j - 2); ++k) {
					mc_blas_vector_at(y, iy) = mc_caddl(mc_blas_vector_at(y, iy), mc_cmull(temp1, mc_blas_vector_at(ap, k)));
					temp2                    = mc_caddl(temp2, mc_cmull(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, ix)));
					ix                       = ix + incx;
					iy                       = iy + incy;
				}
				mc_blas_vector_at(y, jy) = mc_caddl(mc_blas_vector_at(y, jy), mc_caddl(mc_cmull(temp1, mc_blas_vector_at(ap, kk + j - 1)), mc_cmull(alpha, temp2)));
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + j;
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
				mc_blas_vector_at(y, j) = mc_caddl(mc_blas_vector_at(y, j), mc_cmull(temp1, mc_blas_vector_at(ap, kk)));
				k                       = kk + 1;
				for (i = j + 1; i <= n; ++i) {
					mc_blas_vector_at(y, i) = mc_caddl(mc_blas_vector_at(y, i), mc_cmull(temp1, mc_blas_vector_at(ap, k)));
					temp2                   = mc_caddl(temp2, mc_cmull(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, i)));
					k                       = k + 1;
				}
				mc_blas_vector_at(y, j) = mc_caddl(mc_blas_vector_at(y, j), mc_cmull(alpha, temp2));
				kk                      = kk + (n - j + 1);
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
				mc_blas_vector_at(y, jy) = mc_caddl(mc_blas_vector_at(y, jy), mc_cmull(temp1, mc_blas_vector_at(ap, kk)));
				ix                       = jx;
				iy                       = jy;
				for (k = (kk + 1); k <= (kk + n - j); ++k) {
					ix                       = ix + incx;
					iy                       = iy + incy;
					mc_blas_vector_at(y, iy) = mc_caddl(mc_blas_vector_at(y, iy), mc_cmull(temp1, mc_blas_vector_at(ap, k)));
					temp2                    = mc_caddl(temp2, mc_cmull(mc_blas_vector_at(ap, k), mc_blas_vector_at(x, ix)));
				}
				mc_blas_vector_at(y, jy) = mc_caddl(mc_blas_vector_at(y, jy), mc_cmull(alpha, temp2));
				jx                       = jx + incx;
				jy                       = jy + incy;
				kk                       = kk + (n - j + 1);
			}
		}
	}
}

#endif /* !MC_BLAS_SPMV_H */

/* EOF */