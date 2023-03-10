//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_rot.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?rot applies a plane rotation.
 *
 * \synopsis
 *    void ?rot(n, x, incx, y, incy, c, s)
 *    real-floating c, s
 *    int           incx, incy, n
 *    real-floating x(*), y(*)
 *
 * \purpose
 *    ?rot applies a plane rotation matrix to a real sequence of ordered pairs.
 *    If coefficients c and s satisfy c+s=1, the rotation matrix is orthogonal,
 *    and the transformation is called a Givens plane rotation.
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in the input vector `x` and `y`.
 *
 *    [out] x    - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [out] y    - real-floating array of size at least (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [in]  c    - real-floating. Specifies the cosine of the angle of rotation (Givens rotation matrix).
 *    [in]  s    - real-floating. Specifies the sine of the angle of rotation (Givens rotation matrix).
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Linpack.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/details/math/mc_cadd.h>
#include <macadam/details/math/mc_cmul.h>
#include <macadam/details/math/mc_csub.h>

#ifndef MC_BLAS_ROT_H
#define MC_BLAS_ROT_H

#pragma mark - mc_blas_srot -

MC_TARGET_FUNC void mc_blas_srot(const int n, float * x, const int incx, float * y, const int incy, const float c, const float s)
{
	float temp;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = c * mc_blas_vector_at(x, i) + s * mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i) = c * mc_blas_vector_at(y, i) - s * mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = c * mc_blas_vector_at(x, ix) + s * mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = c * mc_blas_vector_at(y, iy) - s * mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_drot -

MC_TARGET_FUNC void mc_blas_drot(const int n, double * x, const int incx, double * y, const int incy, const double c, const double s)
{
	double temp;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = c * mc_blas_vector_at(x, i) + s * mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i) = c * mc_blas_vector_at(y, i) - s * mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = c * mc_blas_vector_at(x, ix) + s * mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = c * mc_blas_vector_at(y, iy) - s * mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_lrot -

MC_TARGET_FUNC void mc_blas_lrot(const int n, long double * x, const int incx, long double * y, const int incy, const long double c, const long double s)
{
	long double temp;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}
	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = c * mc_blas_vector_at(x, i) + s * mc_blas_vector_at(y, i);
			mc_blas_vector_at(y, i) = c * mc_blas_vector_at(y, i) - s * mc_blas_vector_at(x, i);
			mc_blas_vector_at(x, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = c * mc_blas_vector_at(x, ix) + s * mc_blas_vector_at(y, iy);
			mc_blas_vector_at(y, iy) = c * mc_blas_vector_at(y, iy) - s * mc_blas_vector_at(x, ix);
			mc_blas_vector_at(x, ix) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

/* \name
 *    ?rot applies a plane rotation.
 *
 * \synopsis
 *    void ?rot(n, x, incx, y, incy, c, s)
 *    real-floating c, s
 *    int           incx, incy, n
 *    complex       x(*), y(*)
 *
 * \purpose
 *    ?rot applies a plane rotation matrix to a real sequence of ordered pairs.
 *    If coefficients c and s satisfy c+s=1, the rotation matrix is orthogonal,
 *    and the transformation is called a Givens plane rotation.
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in the input vector `x` and `y`.
 *
 *    [out] x    - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [out] y    - complex array of size at least (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [in]  c    - real-floating. Specifies the cosine of the angle of rotation (Givens rotation matrix).
 *    [in]  s    - real-floating. Specifies the sine of the angle of rotation (Givens rotation matrix).
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Linpack.
 */

#pragma mark - mc_blas_csrot -

MC_TARGET_FUNC void mc_blas_csrot(const int n, mc_complex_float_t * x, const int incx, mc_complex_float_t * y, const int incy, const float c, const float s)
{
	mc_complex_float_t temp, cc, cs;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}

	cc = mc_cmplxf(c, 0.0f);
	cs = mc_cmplxf(s, 0.0f);

	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = mc_caddf(mc_cmulf(cc, mc_blas_vector_at(x, i)), mc_cmulf(cs, mc_blas_vector_at(y, i)));
			mc_blas_vector_at(y, i) = mc_csubf(mc_cmulf(cc, mc_blas_vector_at(y, i)), mc_cmulf(cs, mc_blas_vector_at(x, i)));
			mc_blas_vector_at(x, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_caddf(mc_cmulf(cc, mc_blas_vector_at(x, ix)), mc_cmulf(cs, mc_blas_vector_at(y, iy)));
			mc_blas_vector_at(y, iy) = mc_csubf(mc_cmulf(cc, mc_blas_vector_at(y, iy)), mc_cmulf(cs, mc_blas_vector_at(x, ix)));
			mc_blas_vector_at(x, ix) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_zdrot -

MC_TARGET_FUNC void mc_blas_zdrot(const int n, mc_complex_double_t * x, const int incx, mc_complex_double_t * y, const int incy, const double c, const double s)
{
	mc_complex_double_t temp, zc, zs;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}

	zc = mc_cmplx(c, 0.0);
	zs = mc_cmplx(s, 0.0);

	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = mc_cadd(mc_cmul(zc, mc_blas_vector_at(x, i)), mc_cmul(zs, mc_blas_vector_at(y, i)));
			mc_blas_vector_at(y, i) = mc_csub(mc_cmul(zc, mc_blas_vector_at(y, i)), mc_cmul(zs, mc_blas_vector_at(x, i)));
			mc_blas_vector_at(x, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_cadd(mc_cmul(zc, mc_blas_vector_at(x, ix)), mc_cmul(zs, mc_blas_vector_at(y, iy)));
			mc_blas_vector_at(y, iy) = mc_csub(mc_cmul(zc, mc_blas_vector_at(y, iy)), mc_cmul(zs, mc_blas_vector_at(x, ix)));
			mc_blas_vector_at(x, ix) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#pragma mark - mc_blas_qlrot -

MC_TARGET_FUNC void mc_blas_qlrot(const int n, mc_complex_long_double_t * x, const int incx, mc_complex_long_double_t * y, const int incy, const long double c, const long double s)
{
	mc_complex_long_double_t temp, qc, qs;
	int i, ix, iy;

	if (n <= 0) {
		return;
	}

	qc = mc_cmplxl(c, 0.0L);
	qs = mc_cmplxl(s, 0.0L);

	if (incx == 1 && incy == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                    = mc_caddl(mc_cmull(qc, mc_blas_vector_at(x, i)), mc_cmull(qs, mc_blas_vector_at(y, i)));
			mc_blas_vector_at(y, i) = mc_csubl(mc_cmull(qc, mc_blas_vector_at(y, i)), mc_cmull(qs, mc_blas_vector_at(x, i)));
			mc_blas_vector_at(x, i) = temp;
		}
	} else {
		ix = 1;
		iy = 1;
		if (incx < 0) {
			ix = (-(n) + 1) * incx + 1;
		}
		if (incy < 0) {
			iy = (-(n) + 1) * incy + 1;
		}
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 1; i <= n; ++i) {
			temp                     = mc_caddl(mc_cmull(qc, mc_blas_vector_at(x, ix)), mc_cmull(qs, mc_blas_vector_at(y, iy)));
			mc_blas_vector_at(y, iy) = mc_csubl(mc_cmull(qc, mc_blas_vector_at(y, iy)), mc_cmull(qs, mc_blas_vector_at(x, ix)));
			mc_blas_vector_at(x, ix) = temp;
			ix                       = ix + incx;
			iy                       = iy + incy;
		}
	}
}

#endif /* !MC_BLAS_ROT_H */

/* EOF */