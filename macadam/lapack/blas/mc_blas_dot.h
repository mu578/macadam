//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_dot.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?dot computes the inner product of two vectors.
 *
 * \synopsis
 *    real-floating ?dot(n, x, incx, y, incy)
 *    int           incx, incy, n
 *    real-floating x(*), y(*)
 *
 * \purpose
 *    ?dot computes the inner product of two vectors.
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vectors x and y.
 *    [in] x     - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *    [in] y     - real-floating arrays of size at least (1+(n-1)*abs(incy)).
 *    [in] incy  - int. Specifies the increment for the elements of `y`. incy must not be zero.
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
#include <macadam/details/math/mc_conj.h>
#include <macadam/details/math/mc_cmul.h>

#ifndef MC_BLAS_DOT_H
#define MC_BLAS_DOT_H

#pragma mark - mc_blas_sdot -

MC_TARGET_FUNC float mc_blas_sdot(const int n, const float * x, const int incx, const float * y , const int incy)
{
	int i, ix, iy, m, mp1;
	float temp;

	temp = 0.0f;
	if (n <= 0) {
		return temp;
	}
	if (incx == 1 && incy == 1) {
		m = n % 5;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				temp = temp + (mc_blas_vector_at(x, i) * mc_blas_vector_at(y, i));
			}
			if (n < 5) {
				return temp;
			}
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 5) {
			temp = temp + (
				  (mc_blas_vector_at(x, i    ) * mc_blas_vector_at(y, i    ))
				+ (mc_blas_vector_at(x, i + 1) * mc_blas_vector_at(y, i + 1))
				+ (mc_blas_vector_at(x, i + 2) * mc_blas_vector_at(y, i + 2))
				+ (mc_blas_vector_at(x, i + 3) * mc_blas_vector_at(y, i + 3))
				+ (mc_blas_vector_at(x, i + 4) * mc_blas_vector_at(y, i + 4))
			);
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
			temp = temp + (mc_blas_vector_at(x, ix) * mc_blas_vector_at(y, iy));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#pragma mark - mc_blas_dsdot -

MC_TARGET_FUNC double mc_blas_dsdot(const int n, const float * x, const int incx, const float * y , const int incy)
{
	int i, ix, iy, m, mp1;
	double temp;

	temp = 0.0;
	if (n <= 0) {
		return temp;
	}
	if (incx == 1 && incy == 1) {
		m = n % 5;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				temp = temp + (mc_cast(double, mc_blas_vector_at(x, i)) * mc_cast(double, mc_blas_vector_at(y, i)));
			}
			if (n < 5) {
				return temp;
			}
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 5) {
			temp = temp + (
				  (mc_cast(double, mc_blas_vector_at(x, i    )) * mc_cast(double, mc_blas_vector_at(y, i    )))
				+ (mc_cast(double, mc_blas_vector_at(x, i + 1)) * mc_cast(double, mc_blas_vector_at(y, i + 1)))
				+ (mc_cast(double, mc_blas_vector_at(x, i + 2)) * mc_cast(double, mc_blas_vector_at(y, i + 2)))
				+ (mc_cast(double, mc_blas_vector_at(x, i + 3)) * mc_cast(double, mc_blas_vector_at(y, i + 3)))
				+ (mc_cast(double, mc_blas_vector_at(x, i + 4)) * mc_cast(double, mc_blas_vector_at(y, i + 4)))
			);
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
			temp = temp + (mc_cast(double, mc_blas_vector_at(x, ix)) * mc_cast(double, mc_blas_vector_at(y, iy)));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#pragma mark - mc_blas_sdsdot -

MC_TARGET_FUNC float mc_blas_sdsdot(const int n, float b, const float * x, const int incx, const float * y , const int incy)
{
	return mc_cast(float, (mc_cast(double, b) + mc_blas_dsdot(n, x, incx, y, incy)));
}

#pragma mark - mc_blas_ddot -

MC_TARGET_FUNC double mc_blas_ddot(const int n, const double * x, const int incx, const double * y , const int incy)
{
	int i, ix, iy, m, mp1;
	double temp;

	temp = 0.0;
	if (n <= 0) {
		return temp;
	}
	if (incx == 1 && incy == 1) {
		m = n % 5;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				temp = temp + (mc_blas_vector_at(x, i) * mc_blas_vector_at(y, i));
			}
			if (n < 5) {
				return temp;
			}
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 5) {
			temp = temp + (
				  (mc_blas_vector_at(x, i    ) * mc_blas_vector_at(y, i    ))
				+ (mc_blas_vector_at(x, i + 1) * mc_blas_vector_at(y, i + 1))
				+ (mc_blas_vector_at(x, i + 2) * mc_blas_vector_at(y, i + 2))
				+ (mc_blas_vector_at(x, i + 3) * mc_blas_vector_at(y, i + 3))
				+ (mc_blas_vector_at(x, i + 4) * mc_blas_vector_at(y, i + 4))
			);
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
			temp = temp + (mc_blas_vector_at(x, ix) * mc_blas_vector_at(y, iy));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#pragma mark - mc_blas_ldot -

MC_TARGET_FUNC long double mc_blas_ldot(const int n, const long double * x, const int incx, const long double * y , const int incy)
{
	int i, ix, iy, m, mp1;
	long double temp;

	temp = 0.0L;
	if (n <= 0) {
		return temp;
	}
	if (incx == 1 && incy == 1) {
		m = n % 5;
		if (m != 0) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
			for (i = 1; i <= m; ++i) {
				temp = temp + (mc_blas_vector_at(x, i) * mc_blas_vector_at(y, i));
			}
			if (n < 5) {
				return temp;
			}
		}
		mp1 = m + 1;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = mp1; i <= n; i += 5) {
			temp = temp + (
				  (mc_blas_vector_at(x, i    ) * mc_blas_vector_at(y, i    ))
				+ (mc_blas_vector_at(x, i + 1) * mc_blas_vector_at(y, i + 1))
				+ (mc_blas_vector_at(x, i + 2) * mc_blas_vector_at(y, i + 2))
				+ (mc_blas_vector_at(x, i + 3) * mc_blas_vector_at(y, i + 3))
				+ (mc_blas_vector_at(x, i + 4) * mc_blas_vector_at(y, i + 4))
			);
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
			temp = temp + (mc_blas_vector_at(x, ix) * mc_blas_vector_at(y, iy));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

/* \name
 *    ?dotc computes the inner product of two vectors: ?dotc=x_*y.
 *
 * \synopsis
 *    complex ?dotc(n, x, incx, y, incy)
 *    int     incx, incy, n
 *    complex x(*), y(*)
 *
 * \purpose
 *    ?dot computes the inner product of two vectors: ?dotc=x_*y.
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vectors x and y.
 *    [in] x     - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *    [in] y     - complex array of size at least (1+(n-1)*abs(incy)).
 *    [in] incy  - int. Specifies the increment for the elements of `y`. incy must not be zero.
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

#pragma mark - mc_blas_cdotc -

MC_TARGET_FUNC mc_complex_float_t mc_blas_cdotc(const int n, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y , const int incy)
{
	int i, ix, iy;
	mc_complex_float_t temp;

	temp = mc_cmplxf(0.0f, 0.0f);
	if (n <= 0) {
		return temp;
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
			temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_vector_at(x, i)), mc_blas_vector_at(y, i)));
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
			temp = mc_caddf(temp, mc_cmulf(mc_conjf(mc_blas_vector_at(x, ix)), mc_blas_vector_at(y, iy)));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#pragma mark - mc_blas_zdotc -

MC_TARGET_FUNC mc_complex_double_t mc_blas_zdotc(const int n, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y , const int incy)
{
	int i, ix, iy;
	mc_complex_double_t temp;

	temp = mc_cmplx(0.0, 0.0);
	if (n <= 0) {
		return temp;
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
			temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_vector_at(x, i)), mc_blas_vector_at(y, i)));
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
			temp = mc_cadd(temp, mc_cmul(mc_conj(mc_blas_vector_at(x, ix)), mc_blas_vector_at(y, iy)));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#pragma mark - mc_blas_qdotc -

MC_TARGET_FUNC mc_complex_long_double_t mc_blas_qdotc(const int n, const mc_complex_long_double_t * x, const int incx, const mc_complex_long_double_t * y , const int incy)
{
	int i, ix, iy;
	mc_complex_long_double_t temp;

	temp = mc_cmplxl(0.0L, 0.0L);
	if (n <= 0) {
		return temp;
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
			temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_vector_at(x, i)), mc_blas_vector_at(y, i)));
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
			temp = mc_caddl(temp, mc_cmull(mc_conjl(mc_blas_vector_at(x, ix)), mc_blas_vector_at(y, iy)));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

/* \name
 *    ?dotu computes the inner product of two vectors: ?dotu=x'*y.
 *
 * \synopsis
 *    complex ?dotu(n, x, incx, y, incy)
 *    int     incx, incy, n
 *    complex x(*), y(*)
 *
 * \purpose
 *    ?dot computes the inner product of two vectors: ?dotu=x'*y.
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vectors x and y.
 *    [in] x     - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *    [in] y     - complex array of size at least (1+(n-1)*abs(incy)).
 *    [in] incy  - int. Specifies the increment for the elements of `y`. incy must not be zero.
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

#pragma mark - mc_blas_cdotu -

MC_TARGET_FUNC mc_complex_float_t mc_blas_cdotu(const int n, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y , const int incy)
{
	int i, ix, iy;
	mc_complex_float_t temp;

	temp = mc_cmplxf(0.0f, 0.0f);
	if (n <= 0) {
		return temp;
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
			temp = mc_caddf(temp, mc_cmulf(mc_blas_vector_at(x, i), mc_blas_vector_at(y, i)));
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
			temp = mc_caddf(temp, mc_cmulf(mc_blas_vector_at(x, ix), mc_blas_vector_at(y, iy)));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#pragma mark - mc_blas_zdotu -

MC_TARGET_FUNC mc_complex_double_t mc_blas_zdotu(const int n, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y , const int incy)
{
	int i, ix, iy;
	mc_complex_double_t temp;

	temp = mc_cmplx(0.0, 0.0);
	if (n <= 0) {
		return temp;
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
			temp = mc_cadd(temp, mc_cmul(mc_blas_vector_at(x, i), mc_blas_vector_at(y, i)));
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
			temp = mc_cadd(temp, mc_cmul(mc_blas_vector_at(x, ix), mc_blas_vector_at(y, iy)));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#pragma mark - mc_blas_qdotu -

MC_TARGET_FUNC mc_complex_long_double_t mc_blas_qdotu(const int n, const mc_complex_long_double_t * x, const int incx, const mc_complex_long_double_t * y , const int incy)
{
	int i, ix, iy;
	mc_complex_long_double_t temp;

	temp = mc_cmplxl(0.0L, 0.0L);
	if (n <= 0) {
		return temp;
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
			temp = mc_caddl(temp, mc_cmull(mc_blas_vector_at(x, i), mc_blas_vector_at(y, i)));
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
			temp = mc_caddl(temp, mc_cmull(mc_blas_vector_at(x, ix), mc_blas_vector_at(y, iy)));
			ix   = ix + incx;
			iy   = iy + incy;
		}
	}
	return temp;
}

#endif /* !MC_BLAS_DOT_H */

/* EOF */