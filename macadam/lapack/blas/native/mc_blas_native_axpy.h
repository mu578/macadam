//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_axpy.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \synopsis
 *    void ?axpy(n, a, x, incx, y, incy)
 *    int           incx, incy, n
 *    real-floating a, x(*), y(*)
 *
 * \purpose
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in the input vectors x and y.
 *    [in]  a    - real-floating. Specifies the scalar alpha.
 *    [in]  x    - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *    [out] y    - real-floating arrays of size at least (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the increment for the elements of `y`. incy must not be zero.
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

#include <macadam/details/mc_target.h>

#ifndef MC_BLAS_NATIVE_AXPY_H
#define MC_BLAS_NATIVE_AXPY_H

#pragma mark - mc_blas_native_saxpy -

MC_TARGET_FUNC void mc_blas_native_saxpy(const int n, const float a, const float * x, const int incx, float * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_saxpy(n, a, x, incx, y, incy);
#	else
	cblas_saxpy(n, a, x, incx, y, incy);
#	endif
}

#pragma mark - mc_blas_native_daxpy -

MC_TARGET_FUNC void mc_blas_native_daxpy(const int n, const double a, const double * x, const int incx, double * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_daxpy(n, a, x, incx, y, incy);
#	else
	cblas_daxpy(n, a, x, incx, y, incy);
#	endif
}

/* \name
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \synopsis
 *    void ?axpy(n, a, x, incx, y, incy)
 *    int     incx, incy, n
 *    complex a, x(*), y(*)
 *
 * \purpose
 *    ?axpy constant times a vector plus a vector: y=y+(a*x).
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in the input vectors x and y.
 *    [in]  a    - complex. Specifies the scalar alpha.
 *    [in]  x    - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *    [out] y    - complex array of size at least (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the increment for the elements of `y`. incy must not be zero.
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

#pragma mark - mc_blas_native_caxpy -

MC_TARGET_FUNC void mc_blas_native_caxpy(const int n, const mc_complex_float_t a, const mc_complex_float_t * x, const int incx, mc_complex_float_t * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_caxpy(n, &a, x, incx, y, incy);
#	else
	cblas_caxpy(n, &a, x, incx, y, incy);
#	endif
}

#pragma mark - mc_blas_native_zaxpy -

MC_TARGET_FUNC void mc_blas_native_zaxpy(const int n, const mc_complex_double_t a, const mc_complex_double_t * x, const int incx, mc_complex_double_t * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_zaxpy(n, &a, x, incx, y, incy);
#	else
	cblas_zaxpy(n, &a, x, incx, y, incy);
#	endif
}

#endif /* !MC_BLAS_NATIVE_AXPY_H */

/* EOF */