//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_rotg.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?rotg computes the parameters for a Givens rotation.
 *
 * \synopsis
 *    void ?rotg(a, b, c, s)
 *    real-floating a, b, c, s
 *
 * \purpose
 *    ?rotg computes the parameters for a Givens rotation. Given the Cartesian coordinates (a, b)
 *    of a point, returns the parameters c, s, r, and z associated with the Givens rotation.
 *
 * \parameters
 *    [out] a - real-floating. Provides the x-coordinate of the point p, a is overwritten by
 *    the parameter r associated with the Givens rotation.
 *
 *    [out] b - real-floating. Provides the y-coordinate of the point p, b is overwritten by
 *    the parameter z associated with the Givens rotation.
 *
 *    [out] c - real-floating. Contains the cosine associated with the Givens rotation.
 *    [out] s - real-floating. Contains the sine associated with the Givens rotation.
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

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_cabs.h>
#include <macadam/details/math/mc_cdiv.h>
#include <macadam/details/math/mc_cmul.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_BLAS_ROTG_H
#define MC_BLAS_ROTG_H

#pragma mark - mc_blas_srotg -

MC_TARGET_FUNC void mc_blas_srotg(float * a, float * b, float * c, float * s)
{
	float r, roe, scale, z;

	roe = *b;
	if (mc_fabsf(*a) > mc_fabsf(*b)) {
		roe = *a;
	}
	scale = mc_fabsf(*a) + mc_fabsf(*b);
	if (scale == 0.0f) {
		(*c) = 1.0f;
		(*s) = 0.0f;
		r    = 0.0f;
		z    = 0.0f;
	} else {
		r    = scale * mc_sqrtf(mc_raise2f(*a / scale) + mc_raise2f(*b / scale));
		r    = mc_copysignf(1.0f, roe) * r;
		(*c) = *a / r;
		(*s) = *b / r;
		z    = 1.0f;
		if (mc_fabsf(*a) > mc_fabsf(*b)) {
			z = *s;
		}
		if (mc_fabsf(*b) >= mc_fabsf(*a) && *c != 0.0f) {
			z = 1.0f / *c;
		}
	}
	*a = r;
	*b = z;
}

#pragma mark - mc_blas_drotg -

MC_TARGET_FUNC void mc_blas_drotg(double * a, double * b, double * c, double * s)
{
	double r, roe, scale, z;

	roe = *b;
	if (mc_fabs(*a) > mc_fabs(*b)) {
		roe = *a;
	}
	scale = mc_fabs(*a) + mc_fabs(*b);
	if (scale == 0.0) {
		(*c) = 1.0;
		(*s) = 0.0;
		r    = 0.0;
		z    = 0.0;
	} else {
		r    = scale * mc_sqrt(mc_raise2(*a / scale) + mc_raise2(*b / scale));
		r    = mc_copysign(1.0, roe) * r;
		(*c) = *a / r;
		(*s) = *b / r;
		z    = 1.0;
		if (mc_fabs(*a) > mc_fabs(*b)) {
			z = *s;
		}
		if (mc_fabs(*b) >= mc_fabs(*a) && *c != 0.0) {
			z = 1.0 / *c;
		}
	}
	*a = r;
	*b = z;
}

#pragma mark - mc_blas_lrotg -

MC_TARGET_FUNC void mc_blas_lrotg(long double * a, long double * b, long double * c, long double * s)
{
	long double r, roe, scale, z;

	roe = *b;
	if (mc_fabsl(*a) > mc_fabsl(*b)) {
		roe = *a;
	}
	scale = mc_fabsl(*a) + mc_fabsl(*b);
	if (scale == 0.0L) {
		(*c) = 1.0L;
		(*s) = 0.0L;
		r    = 0.0L;
		z    = 0.0L;
	} else {
		r    = scale * mc_sqrtl(mc_raise2l(*a / scale) + mc_raise2l(*b / scale));
		r    = mc_copysignl(1.0L, roe) * r;
		(*c) = *a / r;
		(*s) = *b / r;
		z    = 1.0L;
		if (mc_fabsl(*a) > mc_fabsl(*b)) {
			z = *s;
		}
		if (mc_fabsl(*b) >= mc_fabsl(*a) && *c != 0.0L) {
			z = 1.0L / *c;
		}
	}
	*a = r;
	*b = z;
}

/* \name
 *    ?rotg computes the parameters for a Givens rotation.
 *
 * \synopsis
 *    void ?rotg(a, b, c, s)
 *    complex       a, b, s
 *    real-floating c
 *
 * \purpose
 *    ?rotg computes the parameters for a Givens rotation. Given the Cartesian coordinates (a, b)
 *    of a point, returns the parameters c, s, r, and z associated with the Givens rotation.
 *
 * \parameters
 *    [out] a - complex. Provides the x-coordinate of the point p, a is overwritten by
 *    the parameter r associated with the Givens rotation.
 *
 *    [out] b - complex. Provides the y-coordinate of the point p, b is overwritten by
 *    the parameter z associated with the Givens rotation.
 *
 *    [out] c - real-floating. Contains the cosine associated with the Givens rotation.
 *    [out] s - complex. Contains the sine associated with the Givens rotation.
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

#pragma mark - mc_blas_crotg -

MC_TARGET_FUNC void mc_blas_crotg(mc_complex_float_t * a, mc_complex_float_t * b, float * c, mc_complex_float_t * s)
{
	mc_complex_float_t alpha;
	float norm, scale;

	if (mc_cabsf(*a) == 0.0f) {
		(*c) = 0.0f;
		(*s) = mc_cmplxf(1.0f, 0.0f);
		(*a) = (*b);
	} else {
		scale = mc_cabsf(*a) + mc_cabsf(*b);
		norm  = scale * mc_sqrtf(mc_raise2f(mc_cabsf(mc_cdivf(*a, mc_cmplxf(scale, 0.0f))) + mc_raise2f(mc_cabsf(mc_cdivf(*b, mc_cmplxf(scale, 0.0f))))));
		alpha = mc_cdivf(*a, mc_cmplxf(mc_cabsf(*a), 0.0f));
		(*c)  = mc_cabsf(*a) / norm;
		(*s)  = mc_cdivf(mc_cmulf(alpha, mc_conjf(*b)), mc_cmplxf(norm, 0.0f));
		(*a)  = mc_cmulf(alpha, mc_cmplxf(norm, 0.0f));
	}
}

#pragma mark - mc_blas_zrotg -

MC_TARGET_FUNC void mc_blas_zrotg(mc_complex_double_t * a, mc_complex_double_t * b, double * c, mc_complex_double_t * s)
{
	mc_complex_double_t alpha;
	double norm, scale;

	if (mc_cabs(*a) == 0.0) {
		(*c) = 0.0;
		(*s) = mc_cmplx(1.0, 0.0);
		(*a) = (*b);
	} else {
		scale = mc_cabs(*a) + mc_cabs(*b);
		norm  = scale * mc_sqrt(mc_raise2(mc_cabs(mc_cdiv(*a, mc_cmplx(scale, 0.0))) + mc_raise2(mc_cabs(mc_cdiv(*b, mc_cmplx(scale, 0.0))))));
		alpha = mc_cdiv(*a, mc_cmplx(mc_cabs(*a), 0.0));
		(*c)  = mc_cabs(*a) / norm;
		(*s)  = mc_cdiv(mc_cmul(alpha, mc_conj(*b)), mc_cmplx(norm, 0.0));
		(*a)  = mc_cmul(alpha, mc_cmplx(norm, 0.0));
	}
}

#pragma mark - mc_blas_qrotg -

MC_TARGET_FUNC void mc_blas_qrotg(mc_complex_long_double_t * a, mc_complex_long_double_t * b, long double * c, mc_complex_long_double_t * s)
{
	mc_complex_long_double_t alpha;
	long double norm, scale;

	if (mc_cabsl(*a) == 0.0L) {
		(*c) = 0.0L;
		(*s) = mc_cmplxl(1.0L, 0.0L);
		(*a) = (*b);
	} else {
		scale = mc_cabsl(*a) + mc_cabsl(*b);
		norm  = scale * mc_sqrtl(mc_raise2l(mc_cabsl(mc_cdivl(*a, mc_cmplxl(scale, 0.0L))) + mc_raise2l(mc_cabsl(mc_cdivl(*b, mc_cmplxl(scale, 0.0L))))));
		alpha = mc_cdivl(*a, mc_cmplxl(mc_cabsl(*a), 0.0L));
		(*c)  = mc_cabsl(*a) / norm;
		(*s)  = mc_cdivl(mc_cmull(alpha, mc_conjl(*b)), mc_cmplxl(norm, 0.0L));
		(*a)  = mc_cmull(alpha, mc_cmplxl(norm, 0.0L));
	}
}

#endif /* !MC_BLAS_ROTG_H */

/* EOF */