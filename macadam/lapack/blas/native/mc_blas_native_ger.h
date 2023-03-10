//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_sger.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?ger performs the rank 1 operation: a=alpha*x*y' + a.
 *
 * \synopsis
 *    void ?ger(m, n, alpha, x, incx, y, incy, a, lda)
 *    real-floating alpha
 *    int            incx, incy, lda, m, n
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?ger performs the rank 1 operation: a=alpha*x*y' + a where alpha is a scalar,
 *    x is an m element vector, y is an n element vector and a is an m by n matrix.
 *
 * \parameters
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in]  x     - real-floating array of dimension (at least) (1+(m-1)*abs(incx)). The incremented
 *    array `x` must contain the m element vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y     - real-floating array of dimension (at least) (1+(n-1)*abs(incy)). The incremented
 *    array `y` must contain the n element vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [out] a     - real-floating array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients. a is overwritten by the updated matrix.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
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

#include <macadam/details/mc_target.h>

#ifndef MC_BLAS_NATIVE_GER_H
#define MC_BLAS_NATIVE_GER_H

#pragma mark - mc_blas_native_sger -

MC_TARGET_FUNC void mc_blas_native_sger(const int m, const int n, const float alpha, const float * x, const int incx, const float * y, const int incy, float * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

#	if MC_TARGET_CPP98
	::cblas_sger(ord, m, n, alpha, x, incx, y, incy, a, lda);
#	else
	cblas_sger(ord, m, n, alpha, x, incx, y, incy, a, lda);
#	endif
}

#pragma mark - mc_blas_native_dger -

MC_TARGET_FUNC void mc_blas_native_dger(const int m, const int n, const double alpha, const double * x, const int incx, const double * y, const int incy, double * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

#	if MC_TARGET_CPP98
	::cblas_dger(ord, m, n, alpha, x, incx, y, incy, a, lda);
#	else
	cblas_dger(ord, m, n, alpha, x, incx, y, incy, a, lda);
#	endif
}

/* \name
 *    ?gerc performs the rank 1 operation: a=alpha*x*y_ + a.
 *
 * \synopsis
 *    void ?gerc(m, n, alpha, x, incx, y, incy, a, lda)
 *    complex alpha
 *    int     incx, incy, lda, m, n
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?ger performs the rank 1 operation: a=alpha*x*y_ + a where alpha is a scalar,
 *    x is an m element vector, y is an n element vector and a is an m by n matrix.
 *
 * \parameters
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - complex. Specifies the scalar alpha.
 *
 *    [in]  x     - complex array of dimension (at least) (1+(m-1)*abs(incx)). The incremented
 *    array `x` must contain the m element vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y     - complex array of dimension (at least) (1+(n-1)*abs(incy)). The incremented
 *    array `y` must contain the n element vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [out] a     - complex array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients. a is overwritten by the updated matrix.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
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

#pragma mark - mc_blas_native_cgerc -

MC_TARGET_FUNC void mc_blas_native_cgerc(const int m, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y, const int incy, mc_complex_float_t * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

#	if MC_TARGET_CPP98
	::cblas_cgerc(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	else
	cblas_cgerc(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	endif
}

#pragma mark - mc_blas_native_zgerc -

MC_TARGET_FUNC void mc_blas_native_zgerc(const int m, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y, const int incy, mc_complex_double_t * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

#	if MC_TARGET_CPP98
	::cblas_zgerc(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	else
	cblas_zgerc(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	endif
}

/* \name
 *    ?geru performs the rank 1 operation: a=alpha*x*y' + a.
 *
 * \synopsis
 *    void ?geru(m, n, alpha, x, incx, y, incy, a, lda)
 *    complex alpha
 *    int     incx, incy, lda, m, n
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?ger performs the rank 1 operation: a=alpha*x*y' + a where alpha is a scalar,
 *    x is an m element vector, y is an n element vector and a is an m by n matrix.
 *
 * \parameters
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - complex. Specifies the scalar alpha.
 *
 *    [in]  x     - complex array of dimension (at least) (1+(m-1)*abs(incx)). The incremented
 *    array `x` must contain the m element vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y     - complex array of dimension (at least) (1+(n-1)*abs(incy)). The incremented
 *    array `y` must contain the n element vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [out] a     - complex array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients. a is overwritten by the updated matrix.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
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

#pragma mark - mc_blas_native_cgeru -

MC_TARGET_FUNC void mc_blas_native_cgeru(const int m, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y, const int incy, mc_complex_float_t * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

#	if MC_TARGET_CPP98
	::cblas_cgeru(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	else
	cblas_cgeru(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	endif
}

#pragma mark - mc_blas_native_zgeru -

MC_TARGET_FUNC void mc_blas_native_zgeru(const int m, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y, const int incy, mc_complex_double_t * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

#	if MC_TARGET_CPP98
	::cblas_zgeru(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	else
	cblas_zgeru(ord, m, n, &alpha, x, incx, y, incy, a, lda);
#	endif
}

#endif /* !MC_BLAS_NATIVE_GER_H */

/* EOF */