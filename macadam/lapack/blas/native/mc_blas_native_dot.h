//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_dot.h
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

#include <macadam/details/mc_target.h>

#ifndef MC_BLAS_NATIVE_DOT_H
#define MC_BLAS_NATIVE_DOT_H

#pragma mark - mc_blas_native_sdot -

MC_TARGET_FUNC float mc_blas_native_sdot(const int n, const float * x, const int incx, const float * y , const int incy)
{
#	if MC_TARGET_CPP98
	return ::cblas_sdot(n, x, incx, y, incy);
#	else
	return cblas_sdot(n, x, incx, y, incy);
#	endif
}

#pragma mark - mc_blas_native_ddot -

MC_TARGET_FUNC double mc_blas_native_ddot(const int n, const double * x, const int incx, const double * y , const int incy)
{
#	if MC_TARGET_CPP98
	return ::cblas_ddot(n, x, incx, y, incy);
#	else
	return cblas_ddot(n, x, incx, y, incy);
#	endif
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

MC_TARGET_FUNC mc_complex_float_t mc_blas_native_cdotc(const int n, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y , const int incy)
{
	mc_complex_float_t dotc = mc_cmplxf(0.0f, 0.0f);
#	if MC_TARGET_CPP98
	::cblas_cdotc_sub(n, x, incx, y, incy, &dotc);
#	else
	cblas_cdotc_sub(n, x, incx, y, incy, &dotc);
#	endif
	return dotc;
}

#pragma mark - mc_blas_zdotc -

MC_TARGET_FUNC mc_complex_double_t mc_blas_native_zdotc(const int n, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y , const int incy)
{
	mc_complex_double_t dotc = mc_cmplx(0.0, 0.0);
#	if MC_TARGET_CPP98
	::cblas_zdotc_sub(n, x, incx, y, incy, &dotc);
#	else
	cblas_zdotc_sub(n, x, incx, y, incy, &dotc);
#	endif
	return dotc;
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

MC_TARGET_FUNC mc_complex_float_t mc_blas_native_cdotu(const int n, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y , const int incy)
{
	mc_complex_float_t dotu = mc_cmplxf(0.0f, 0.0f);
#	if MC_TARGET_CPP98
	::cblas_cdotu_sub(n, x, incx, y, incy, &dotu);
#	else
	cblas_cdotu_sub(n, x, incx, y, incy, &dotu);
#	endif
	return dotu;
}

#pragma mark - mc_blas_zdotu -

MC_TARGET_FUNC mc_complex_double_t mc_blas_native_zdotu(const int n, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y , const int incy)
{
	mc_complex_double_t dotu = mc_cmplx(0.0, 0.0);
#	if MC_TARGET_CPP98
	::cblas_zdotu_sub(n, x, incx, y, incy, &dotu);
#	else
	cblas_zdotu_sub(n, x, incx, y, incy, &dotu);
#	endif
	return dotu;
}

#endif /* !MC_BLAS_NATIVE_DOT_H */

/* EOF */