//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_spr2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?spr2 performs the symmetric rank 2 operation:
 *    a=alpha*x*y' + alpha*y*x' + a.
 *
 * \synopsis
 *    void ?spr2(uplo, n, alpha, x, incx, y, incy, ap)
 *    real-floating alpha
 *    int            incx, incy, n
 *    char           uplo
 *    real-floating ap(*), x(*), y(*)
 *
 * \purpose
 *    ?spr2 performs the symmetric rank 2 operation: a=alpha*x*y' + alpha*y*x' + a where alpha is a scalar,
 *    `x` and `y` are n element vectors and `a` is an n by n symmetric matrix, supplied in packed form.
 *
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the matrix `a` is supplied
 *    in the packed array `ap` as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` supplied in `ap`.
 *    uplo='L' or 'l', the lower triangular part of `a` supplied in `ap`.
 *
 *    [in] n     - int. Specifies the order of the symmetric matrix `a`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [int] x    - real-floating array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must
 *    contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [int] y    - real-floating array of size at least (1+(n-1)*abs(incy)). The incremented array `y` must
 *    contain the vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [out] ap   - real-floating array of dimension (at least) ((n*(n+1))/2).
 *    With uplo='U' or 'u', the array `ap` must contain the upper triangular part of the symmetric matrix
 *    packed sequentially, column by column, so that ap(1) contains a(1,1), ap(2) and ap(3) contain a(1,2)
 *    and a(2,2) respectively, and so on.
 *
 *    With uplo='L' or 'l', the array `a`p must contain the lower triangular part of the symmetric matrix
 *    packed sequentially, column by column, so that ap(1) contains a(1,1), ap(2) and ap(3) contain a(2,1)
 *    and a(3,1) respectively, and so on.
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

#ifndef MC_BLAS_SPR2_H
#define MC_BLAS_SPR2_H

#pragma mark - mc_blas_sspr2 -

MC_TARGET_FUNC void mc_blas_sspr2(const char uplo, const int n, const float alpha, const float * x, const int incx, const float * y, const int incy, float * ap)
{
	const float zero = 0.0f;

	float temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_x = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_x = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_x, 'U') && !mc_blas_lsame(uplo_x, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("SSPR2 ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_x, 'U')) {
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
					k     = kk;
					for (i = 1; i <= j; ++i) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
						k                        = k + 1;
					}
				}
				kk = kk +  j;
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
					ix = kx;
					iy = ky;
					for (k = kk; k <= (kk + j - 1); ++k) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                       = ix + incx;
						iy                       = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
				kk = kk +  j;
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
					k     = kk;
					for (i = j; i <= n; ++i) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
						k                        = k + 1;
					}
				}
				kk = kk + n - j + 1;
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
					for (k = kk; k <= (kk + n - j); ++k) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                       = ix + incx;
						iy                       = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
				kk = kk + n - j + 1;
			}
		}
	}
}

#pragma mark - mc_blas_dspr2 -

MC_TARGET_FUNC void mc_blas_dspr2(const char uplo, const int n, const double alpha, const double * x, const int incx, const double * y, const int incy, double * ap)
{
	const double zero = 0.0;

	double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_x = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_x = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_x, 'U') && !mc_blas_lsame(uplo_x, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("DSPR2 ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_x, 'U')) {
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
					k     = kk;
					for (i = 1; i <= j; ++i) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
						k                        = k + 1;
					}
				}
				kk = kk +  j;
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
					ix = kx;
					iy = ky;
					for (k = kk; k <= (kk + j - 1); ++k) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                       = ix + incx;
						iy                       = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
				kk = kk +  j;
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
					k     = kk;
					for (i = j; i <= n; ++i) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
						k                        = k + 1;
					}
				}
				kk = kk + n - j + 1;
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
					for (k = kk; k <= (kk + n - j); ++k) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                       = ix + incx;
						iy                       = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
				kk = kk + n - j + 1;
			}
		}
	}
}

#pragma mark - mc_blas_lspr2 -

MC_TARGET_FUNC void mc_blas_lspr2(const char uplo, const int n, const long double alpha, const long double * x, const int incx, const long double * y, const int incy, long double * ap)
{
	const long double zero = 0.0L;

	long double temp1, temp2;
	int i, info, ix, iy, j, jx, jy, k, kk, kx, ky;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	const char uplo_x = mc_blas_lsame(uplo, 'U') ? 'L' : (mc_blas_lsame(uplo, 'L') ? 'U' : 'D');
#	else
	const char uplo_x = uplo;
#	endif

	info = 0;
	if (!mc_blas_lsame(uplo_x, 'U') && !mc_blas_lsame(uplo_x, 'L')) {
		info = 1;
	} else if (n < 0) {
		info = 2;
	} else if (incx == 0) {
		info = 5;
	} else if (incy == 0) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("LSPR2 ", info);
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

	kk = 1;
	if (mc_blas_lsame(uplo_x, 'U')) {
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
					k     = kk;
					for (i = 1; i <= j; ++i) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
						k                        = k + 1;
					}
				}
				kk = kk +  j;
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
					ix = kx;
					iy = ky;
					for (k = kk; k <= (kk + j - 1); ++k) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                       = ix + incx;
						iy                       = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
				kk = kk +  j;
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
					k     = kk;
					for (i = j; i <= n; ++i) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, i) * temp1 + mc_blas_vector_at(y, i) * temp2;
						k                        = k + 1;
					}
				}
				kk = kk + n - j + 1;
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
					for (k = kk; k <= (kk + n - j); ++k) {
						mc_blas_vector_at(ap, k) = mc_blas_vector_at(ap, k) + mc_blas_vector_at(x, ix) * temp1 + mc_blas_vector_at(y, iy) * temp2;
						ix                       = ix + incx;
						iy                       = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
				kk = kk + n - j + 1;
			}
		}
	}
}

#endif /* !MC_BLAS_SPR2_H */

/* EOF */