//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_gemv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?gemv performs one of the matrix-vector operations:
 *    y=alpha*a*x + beta*y or y=alpha*a'*x + beta*y.
 *
 * \synopsis
 *    void ?gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
 *    real-floating alpha, beta
 *    int           incx, incy, lda, m, n
 *    char          trans
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?gemv performs one of the matrix-vector operations: y=alpha*a*x + beta*y or
 *    y=alpha*a'*x + beta*y where alpha and beta are scalars, `x` and `y` arevectors
 *    and a is an m by n matrix.
 *
 * \parameters
 *    [in]  trans - char. Specifies the operation to be performed as follows:
 *    trans= 'N' or 'n' y=alpha*a*x + beta*y.
 *    trans= 'T' or 't' y=alpha*a'*x + beta*y.
 *    trans= 'C' or 'c' y=alpha*a'*x + beta*y.
 *
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in]  a     - real-floating array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
 *
 *    [in]  x     - real-floating array of dimension (at least) (1+(n-1)*abs(incx)) when trans='N' or 'n'
 *    and at least (1+(m-1)*abs(incx)) otherwise. The incremented array `x` must contain the vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  beta  - real-floating. Specifies the scalar beta. When beta is supplied as zero then y need
 *    not be set on input.

 *    [out] y     - real-floating array of dimension (at least) (1+(m-1)*abs(incy)) when trans='N' or 'n'
 *    and at least (1+(n-1)*abs(incy)) otherwise. With beta non-zero, the incremented array `y` must contain
 *    the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
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

#include <macadam/lapack/blas/mc_blas_lsame.h>

#ifndef MC_BLAS_NATIVE_GEMV_H
#define MC_BLAS_NATIVE_GEMV_H

#pragma mark - mc_blas_native_sgemv -

MC_TARGET_FUNC void mc_blas_native_sgemv(const char trans, const int m, const int n, const float alpha, const float * a, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_sgemv(ord, ta, m, n, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_sgemv(ord, ta, m, n, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

#pragma mark - mc_blas_native_dgemv -

MC_TARGET_FUNC void mc_blas_native_dgemv(const char trans, const int m, const int n, const double alpha, const double * a, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_dgemv(ord, ta, m, n, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_dgemv(ord, ta, m, n, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

/* \name
 *    ?gemv performs one of the matrix-vector operations:
 *    y=alpha*a*x + beta*y or y=alpha*a'*x + beta*y or y=alpha*a_*x + beta*y.
 *
 * \synopsis
 *    void ?gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
 *    complex alpha, beta
 *    int     incx, incy, lda, m, n
 *    char    trans
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?gemv performs one of the matrix-vector operations: y=alpha*a*x + beta*y or
 *    y=alpha*a'*x + beta*y or y=alpha*a_*x + beta*y where alpha and beta are scalars,
 *    `x` and `y` arevectors and a is an m by n matrix.
 *
 * \parameters
 *    [in]  trans - char. Specifies the operation to be performed as follows:
 *    trans= 'N' or 'n' y=alpha*a*x + beta*y.
 *    trans= 'T' or 't' y=alpha*a'*x + beta*y.
 *    trans= 'C' or 'c' y=alpha*a_*x + beta*y.
 *
 *    [in]  m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in]  n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *
 *    [in]  alpha - complex. Specifies the scalar alpha.
 *
 *    [in]  a     - complex array of dimension (lda, n), the leading m by n part of the
 *    array a must contain the matrix of coefficients.
 *
 *    [in]  lda   - int. Specifies the first dimension of `a`, lda must be at least max(1, m).
 *
 *    [in]  x     - complex array of dimension (at least) (1+(n-1)*abs(incx)) when trans='N' or 'n'
 *    and at least (1+(m-1)*abs(incx)) otherwise. The incremented array `x` must contain the vector `x`.
 *
 *    [in]  incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  beta  - complex. Specifies the scalar beta. When beta is supplied as zero then y need
 *    not be set on input.

 *    [out] y     - complex array of dimension (at least) (1+(m-1)*abs(incy)) when trans='N' or 'n'
 *    and at least (1+(n-1)*abs(incy)) otherwise. With beta non-zero, the incremented array `y` must contain
 *    the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in]  incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
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

#pragma mark - mc_blas_native_cgemv -

MC_TARGET_FUNC void mc_blas_native_cgemv(const char trans, const int m, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * x, const int incx, const mc_complex_float_t beta, mc_complex_float_t * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_cgemv(ord, ta, m, n, &alpha, a, lda, x, incx, &beta, y, incy);
#	else
	cblas_cgemv(ord, ta, m, n, &alpha, a, lda, x, incx, &beta, y, incy);
#	endif
}

#pragma mark - mc_blas_native_zgemv -

MC_TARGET_FUNC void mc_blas_native_zgemv(const char trans, const int m, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * x, const int incx, const mc_complex_double_t beta, mc_complex_double_t * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_zgemv(ord, ta, m, n, &alpha, a, lda, x, incx, &beta, y, incy);
#	else
	cblas_zgemv(ord, ta, m, n, &alpha, a, lda, x, incx, &beta, y, incy);
#	endif
}

#endif /* !MC_BLAS_NATIVE_GEMV_H */

/* EOF */