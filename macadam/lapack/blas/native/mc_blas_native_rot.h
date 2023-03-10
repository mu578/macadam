//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_rot.h
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

#ifndef MC_BLAS_NATIVE_ROT_H
#define MC_BLAS_NATIVE_ROT_H

#pragma mark - mc_blas_native_srot -

MC_TARGET_FUNC void mc_blas_native_srot(const int n, float * x, const int incx, float * y, const int incy, const float c, const float s)
{
#	if MC_TARGET_CPP98
	::cblas_srot(n, x, incx, y, incy, c, s);
#	else
	cblas_srot(n, x, incx, y, incy, c, s);
#	endif
}

#pragma mark - mc_blas_native_drot -

MC_TARGET_FUNC void mc_blas_native_drot(const int n, double * x, const int incx, double * y, const int incy, const double c, const double s)
{
#	if MC_TARGET_CPP98
	::cblas_drot(n, x, incx, y, incy, c, s);
#	else
	cblas_drot(n, x, incx, y, incy, c, s);
#	endif
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

#pragma mark - mc_blas_native_csrot -

MC_TARGET_FUNC void mc_blas_native_csrot(const int n, mc_complex_float_t * x, const int incx, mc_complex_float_t * y, const int incy, const float c, const float s)
{
#	if !MC_TARGET_BLAS_USE_OPENBLAS
#		if MC_TARGET_CPP98
			::cblas_csrot(n, x, incx, y, incy, c, s);
#		else
			cblas_csrot(n, x, incx, y, incy, c, s);
#		endif
#	else
	mc_unused(n);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(y);
	mc_unused(incy);
	mc_unused(c);
	mc_unused(s);
#	endif
}



#pragma mark - mc_blas_native_zdrot -

MC_TARGET_FUNC void mc_blas_native_zdrot(const int n, mc_complex_double_t * x, const int incx, mc_complex_double_t * y, const int incy, const double c, const double s)
{
#	if !MC_TARGET_BLAS_USE_OPENBLAS
#		if MC_TARGET_CPP98
			::cblas_zdrot(n, x, incx, y, incy, c, s);
#		else
			cblas_zdrot(n, x, incx, y, incy, c, s);
#		endif
#	else
	mc_unused(n);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(y);
	mc_unused(incy);
	mc_unused(c);
	mc_unused(s);
#	endif
}

#endif /* !MC_BLAS_NATIVE_ROT_H */

/* EOF */