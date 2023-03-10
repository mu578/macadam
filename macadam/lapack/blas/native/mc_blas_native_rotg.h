//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_rotg.h
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

#include <macadam/lapack/blas/mc_blas_access.h>

#ifndef MC_BLAS_NATIVE_ROTG_H
#define MC_BLAS_NATIVE_ROTG_H

#pragma mark - mc_blas_native_srotg -

MC_TARGET_FUNC void mc_blas_native_srotg(float * a, float * b, float * c, float * s)
{
#	if MC_TARGET_CPP98
	::cblas_srotg(a, b, c, s);
#	else
	cblas_srotg(a, b, c, s);
#	endif
}

#pragma mark - mc_blas_native_drotg -

MC_TARGET_FUNC void mc_blas_native_drotg(double * a, double * b, double * c, double * s)
{
#	if MC_TARGET_CPP98
	::cblas_drotg(a, b, c, s);
#	else
	cblas_drotg(a, b, c, s);
#	endif
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

#pragma mark - mc_blas_native_crotg -

MC_TARGET_FUNC void mc_blas_native_crotg(mc_complex_float_t * a, mc_complex_float_t * b, float * c, mc_complex_float_t * s)
{
#	if !MC_TARGET_BLAS_USE_OPENBLAS
#		if MC_TARGET_CPP98
			::cblas_crotg(a, b, c, s);
#		else
			cblas_crotg(a, b, c, s);
#		endif
#	else
	mc_unused(a);
	mc_unused(b);
	mc_unused(c);
	mc_unused(s);
#	endif
}

#pragma mark - mc_blas_native_zrotg -

MC_TARGET_FUNC void mc_blas_native_zrotg(mc_complex_double_t * a, mc_complex_double_t * b, double * c, mc_complex_double_t * s)
{
#	if !MC_TARGET_BLAS_USE_OPENBLAS
#		if MC_TARGET_CPP98
			::cblas_zrotg(a, b, c, s);
#		else
			cblas_zrotg(a, b, c, s);
#		endif
#	else
	mc_unused(a);
	mc_unused(b);
	mc_unused(c);
	mc_unused(s);
#	endif
}

#endif /* !MC_BLAS_NATIVE_ROTG_H */

/* EOF */