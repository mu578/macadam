//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_symv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?symv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y.
 *
 * \synopsis
 *    void ?symv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
 *    real-floating alpha, beta
 *    int           incx, incy, lda, n
 *    char          uplo
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?symv performs the matrix-vector operation: c=alpha*a*b + beta*c or c where alpha and
 *    beta are scalars, `x` and `y` are n element vectors and `a` is an n by n symmetric matrix.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `a`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `a` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `a` is to be referenced.
 *
 *    [in] n     - int. Specifies the ord of the matrix `a`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper
 *    triangular part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower
 *    triangular part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least max(1, n).
 *
 *    [in]  x    - real-floating array of dimension (at least) (1+(n-1)*abs(incx)). The incremented array `x`
 *    must contain the vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - real-floating. Specifies the scalar beta. when beta is supplied as zero then `y` need not
 *    be set on input.
 *
 *    [out] y    - real-floating array of dimension (at least) (1+(n-1)*abs(incy)). The incremented array `y`
 *    must contain the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`. incy must not be zero.
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

#ifndef MC_BLAS_NATIVE_SYMV_H
#define MC_BLAS_NATIVE_SYMV_H

#pragma mark - mc_blas_native_ssymv -

MC_TARGET_FUNC void mc_blas_native_ssymv(const char uplo, const int n, const float alpha, const float * a, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_ssymv(ord, ul, n, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_ssymv(ord, ul, n, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

#pragma mark - mc_blas_native_dsymv -

MC_TARGET_FUNC void mc_blas_native_dsymv(const char uplo, const int n, const double alpha, const double * a, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_dsymv(ord, ul, n, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_dsymv(ord, ul, n, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

/* \name
 *    ?symv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y.
 *
 * \synopsis
 *    void ?symv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
 *    complex alpha, beta
 *    int           incx, incy, lda, n
 *    char          uplo
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?symv performs the matrix-vector operation: c=alpha*a*b + beta*c or c where alpha and
 *    beta are scalars, `x` and `y` are n element vectors and `a` is an n by n symmetric matrix.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `a`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `a` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `a` is to be referenced.
 *
 *    [in] n     - int. Specifies the ord of the matrix `a`, n must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper
 *    triangular part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower
 *    triangular part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least max(1, n).
 *
 *    [in]  x    - complex array of dimension (at least) (1+(n-1)*abs(incx)). The incremented array `x`
 *    must contain the vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - complex. Specifies the scalar beta. when beta is supplied as zero then `y` need not
 *    be set on input.
 *
 *    [out] y    - complex array of dimension (at least) (1+(n-1)*abs(incy)). The incremented array `y`
 *    must contain the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`. incy must not be zero.
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

#pragma mark - mc_blas_native_csymv -

MC_TARGET_FUNC void mc_blas_native_csymv(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * x, const int incx, const mc_complex_float_t beta, mc_complex_float_t * y, const int incy)
{
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		if MC_TARGET_BLAS_USE_CLAYOUT
			const enum CBLAS_ORDER ord = CblasRowMajor;
#		else
			const enum CBLAS_ORDER ord = CblasColMajor;
#		endif

		const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#		if MC_TARGET_CPP98
			::cblas_csymv(ord, ul, n, &alpha, a, lda, x, incx, &beta, y, incy);
#		else
			cblas_csymv(ord, ul, n, &alpha, a, lda, x, incx, &beta, y, incy);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(a);
	mc_unused(lda);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(beta);
	mc_unused(y);
	mc_unused(incy);
#	endif
}

#pragma mark - mc_blas_native_zsymv -

MC_TARGET_FUNC void mc_blas_native_zsymv(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * x, const int incx, const mc_complex_double_t beta, mc_complex_double_t * y, const int incy)
{
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		if MC_TARGET_BLAS_USE_CLAYOUT
			const enum CBLAS_ORDER ord = CblasRowMajor;
#		else
			const enum CBLAS_ORDER ord = CblasColMajor;
#		endif

		const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#		if MC_TARGET_CPP98
			::cblas_zsymv(ord, ul, n, &alpha, a, lda, x, incx, &beta, y, incy);
#		else
			cblas_zsymv(ord, ul, n, &alpha, a, lda, x, incx, &beta, y, incy);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(a);
	mc_unused(lda);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(beta);
	mc_unused(y);
	mc_unused(incy);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SYMV_H */

/* EOF */