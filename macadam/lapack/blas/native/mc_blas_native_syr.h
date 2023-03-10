//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_syr.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?syr performs the symmetric rank 1 operation:
 *    a=alpha*x*x' + a.
 *
 * \synopsis
 *    void ?syr(uplo, n, alpha, x, incx, a, lda)
 *    real-floating alpha
 *    int           incx, lda, n
 *    char          uplo
 *    real-floating a(lda, *), x(*)
 *
 * \purpose
 *    ?syr performs the symmetric rank 1 operation: a=alpha*x*x' + a where alpha is a scalar,
 *    `x` is an n element vector and `a` is an n by n symmetric matrix.
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

#ifndef MC_BLAS_NATIVE_SYR_H
#define MC_BLAS_NATIVE_SYR_H

#pragma mark - mc_blas_native_ssyr -

MC_TARGET_FUNC void mc_blas_native_ssyr(const char uplo, const int n, const float alpha, const float * x, const int incx, float * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_ssyr(ord, ul, n, alpha, x, incx, a, lda);
#	else
	cblas_ssyr(ord, ul, n, alpha, x, incx, a, lda);
#	endif
}

#pragma mark - mc_blas_native_dsyr -

MC_TARGET_FUNC void mc_blas_native_dsyr(const char uplo, const int n, const double alpha, const double * x, const int incx, double * a, const int lda)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_dsyr(ord, ul, n, alpha, x, incx, a, lda);
#	else
	cblas_dsyr(ord, ul, n, alpha, x, incx, a, lda);
#	endif
}

/* \name
 *    ?syr performs the symmetric rank 1 operation:
 *    a=alpha*x*x_ + a.
 *
 * \synopsis
 *    void ?syr(uplo, n, alpha, x, incx, a, lda)
 *    complex alpha
 *    int     incx, lda, n
 *    char    uplo
 *    complex a(lda, *), x(*)
 *
 * \purpose
 *    ?syr performs the symmetric rank 1 operation: a=alpha*x*x_ + a where alpha is a scalar,
 *    `x` is an n element vector and `a` is an n by n symmetric matrix.
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

#pragma mark - mc_blas_csyr -

MC_TARGET_FUNC void mc_blas_native_csyr(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, mc_complex_float_t * a, const int lda)
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
			::cblas_csyr(ord, ul, n, &alpha, x, incx, a, lda);
#		else
			cblas_csyr(ord, ul, n, &alpha, x, incx, a, lda);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(a);
	mc_unused(lda);
#	endif
}

#pragma mark - mc_blas_zsyr -

MC_TARGET_FUNC void mc_blas_native_zsyr(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, mc_complex_double_t * a, const int lda)
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
			::cblas_zsyr(ord, ul, n, &alpha, x, incx, a, lda);
#		else
			cblas_zsyr(ord, ul, n, &alpha, x, incx, a, lda);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(a);
	mc_unused(lda);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SYR_H */

/* EOF */