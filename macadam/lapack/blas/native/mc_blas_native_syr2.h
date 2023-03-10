//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_syr2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?syr2 performs the symmetric rank 2 operation:
 *    a=alpha*x*y' + alpha*y*x' + a.
 *
 * \synopsis
 *    void ?syr2(uplo, n, alpha, x, incx, y, incy, a, lda)
 *    real-floating alpha
 *    int           incx, incy, lda, n
 *    char          uplo
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?syr2 performs the symmetric rank 2 operation: a=alpha*x*y' + alpha*y*x' + a where alpha
 *    is a scalar, `x` and `y` are n element vectors and `a` is an n by n symmetric matrix.
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
 *    [in]  x    - real-floating array of dimension (at least) (1+(n-1)*abs(incx)). The incremented array `x`
 *    must contain the vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y    - real-floating array of dimension (at least) (1+(n-1)*abs(incy)). The incremented array `y`
 *    must contain the vector `y`.
 *
 *    [in]  incy - int. Specifies the increment for the elements of `y`, incx must not be zero.
 *
 *    [out] a    - real-floating array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper
 *    triangular part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced.
 *    The upper triangular part of the array `a` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower
 *    triangular part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *    The lower triangular part of the array `a` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least max(1, n).
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

#ifndef MC_BLAS_NATIVE_SYR2_H
#define MC_BLAS_NATIVE_SYR2_H

#pragma mark - mc_blas_native_ssyr2 -

MC_TARGET_FUNC void mc_blas_native_ssyr2(const char uplo, const int n, const float alpha, const float * x, const int incx, const float * y, const int incy, float * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_ssyr2(ord, ul, n, alpha, x, incx, y, incy, a, lda);
#	else
	cblas_ssyr2(ord, ul, n, alpha, x, incx, y, incy, a, lda);
#	endif
}

#pragma mark - mc_blas_native_dsyr2 -

MC_TARGET_FUNC void mc_blas_native_dsyr2(const char uplo, const int n, const double alpha, const double * x, const int incx, const double * y, const int incy, double * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_dsyr2(ord, ul, n, alpha, x, incx, y, incy, a, lda);
#	else
	cblas_dsyr2(ord, ul, n, alpha, x, incx, y, incy, a, lda);
#	endif
}

/* \name
 *    ?syr2 performs the symmetric rank 2 operation:
 *    a=alpha*x*y_ + alpha*y*x_ + a.
 *
 * \synopsis
 *    void ?syr2(uplo, n, alpha, x, incx, y, incy, a, lda)
 *    complex alpha
 *    int     incx, incy, lda, n
 *    char    uplo
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?syr2 performs the symmetric rank 2 operation: a=alpha*x*y' + alpha*y*x' + a where alpha
 *    is a scalar, `x` and `y` are n element vectors and `a` is an n by n symmetric matrix.
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
 *    [in]  x    - complex array of dimension (at least) (1+(n-1)*abs(incx)). The incremented array `x`
 *    must contain the vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in]  y    - complex array of dimension (at least) (1+(n-1)*abs(incy)). The incremented array `y`
 *    must contain the vector `y`.
 *
 *    [in]  incy - int. Specifies the increment for the elements of `y`, incx must not be zero.
 *
 *    [out] a    - complex array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper
 *    triangular part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced.
 *    The upper triangular part of the array `a` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower
 *    triangular part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *    The lower triangular part of the array `a` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least max(1, n).
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

#pragma mark - mc_blas_native_csyr2 -

MC_TARGET_FUNC void mc_blas_native_csyr2(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, const mc_complex_float_t * y, const int incy, mc_complex_float_t * a, const int lda)
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
			::cblas_csyr2(ord, ul, n, &alpha, x, incx, y, incy, a, lda);
#		else
			cblas_csyr2(ord, ul, n, &alpha, x, incx, y, incy, a, lda);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(y);
	mc_unused(incy);
	mc_unused(a);
	mc_unused(lda);
#	endif
}
#pragma mark - mc_blas_native_zsyr2 -

MC_TARGET_FUNC void mc_blas_native_zsyr2(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, const mc_complex_double_t * y, const int incy, mc_complex_double_t * a, const int lda)
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
			::cblas_zsyr2(ord, ul, n, &alpha, x, incx, y, incy, a, lda);
#		else
			cblas_zsyr2(ord, ul, n, &alpha, x, incx, y, incy, a, lda);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(y);
	mc_unused(incy);
	mc_unused(a);
	mc_unused(lda);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SYR2_H */

/* EOF */