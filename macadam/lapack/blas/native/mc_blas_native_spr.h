//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_spr.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?spmv performs the symmetric rank 1 operation:
 *    a=alpha*x*x' + a.
 *
 * \synopsis
 *    void ?spmv(uplo, n, alpha, x, incx, ap)
 *    real-floating alpha
 *    int            incx, n
 *    char           uplo
 *    real-floating ap(*), x(*)
 *
 * \purpose
 *    ?spmv performs the symmetric rank 1 operation: a=alpha*x*x' + a where alpha is a scalar,
 *    `x` is an n element vector and `a` is an n by n symmetric matrix, supplied in packed form.
 *
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the matrix `a` is supplied
 *    in the packed array `ap` as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` supplied in `ap`.
 *    uplo='L' or 'l', the lower triangular part of `a` supplied in `ap`.
 *
 *    [in] n     - int. Specifies the ord of the symmetric matrix `a`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [int] x    - real-floating array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must
 *    contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] ap    - real-floating array of dimension (at least) ((n*(n+1))/2).
 *    With uplo='U' or 'u', the array `ap` must contain the upper triangular part of the symmetric matrix
 *    packed sequentially, column by column, so that ap(1) contains a(1,1), ap(2) and ap(3) contain a(1,2)
 *    and a(2,2) respectively, and so on.
 *
 *    With uplo='L' or 'l', the array `a`p must contain the lower triangular part of the symmetric matrix
 *    packed sequentially, column by column, so that ap(1) contains a(1,1), ap(2) and ap(3) contain a(2,1)
 *    and a(3,1) respectively, and so on.
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

#ifndef MC_BLAS_NATIVE_SPR_H
#define MC_BLAS_NATIVE_SPR_H

#pragma mark - mc_blas_native_sspr -

MC_TARGET_FUNC void mc_blas_native_sspr(const char uplo, const int n, const float alpha, const float * x, const int incx, float * ap)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_sspr(ord, ul, n, alpha, x, incx, ap);
#	else
	cblas_sspr(ord, ul, n, alpha, x, incx, ap);
#	endif
}

#pragma mark - mc_blas_native_dspr -

MC_TARGET_FUNC void mc_blas_native_dspr(const char uplo, const int n, const double alpha, const double * x, const int incx, double * ap)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_dspr(ord, ul, n, alpha, x, incx, ap);
#	else
	cblas_dspr(ord, ul, n, alpha, x, incx, ap);
#	endif
}

/* \name
 *    ?spr performs the symmetric rank 1 operation:
 *    a=alpha*x*x_ + a.
 *
 * \synopsis
 *    void ?spr(uplo, n, alpha, x, incx, ap)
 *    complex alpha
 *    int     incx, n
 *    char    uplo
 *    complex ap(*), x(*)
 *
 * \purpose
 *    ?spr performs the symmetric rank 1 operation: a=alpha*x*x' + a where alpha is a scalar,
 *    `x` is an n element vector and `a` is an n by n symmetric matrix, supplied in packed form.
 *
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the matrix `a` is supplied
 *    in the packed array `ap` as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` supplied in `ap`.
 *    uplo='L' or 'l', the lower triangular part of `a` supplied in `ap`.
 *
 *    [in] n     - int. Specifies the ord of the symmetric matrix `a`, n must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [int] x    - complex array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must
 *    contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [out] ap   - complex array of dimension (at least) ((n*(n+1))/2).
 *    With uplo='U' or 'u', the array `ap` must contain the upper triangular part of the symmetric matrix
 *    packed sequentially, column by column, so that ap(1) contains a(1,1), ap(2) and ap(3) contain a(1,2)
 *    and a(2,2) respectively, and so on.
 *
 *    With uplo='L' or 'l', the array `a`p must contain the lower triangular part of the symmetric matrix
 *    packed sequentially, column by column, so that ap(1) contains a(1,1), ap(2) and ap(3) contain a(2,1)
 *    and a(3,1) respectively, and so on.
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

#pragma mark - mc_blas_native_cspr -

MC_TARGET_FUNC void mc_blas_native_cspr(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * x, const int incx, mc_complex_float_t * ap)
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
			::cblas_cspr(ord, ul, n, &alpha, x, incx, ap);
#		else
			cblas_cspr(ord, ul, n, &alpha, x, incx, ap);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(ap);
#	endif
}

#pragma mark - mc_blas_native_zspr -

MC_TARGET_FUNC void mc_blas_native_zspr(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * x, const int incx, mc_complex_double_t * ap)
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
			::cblas_zspr(ord, ul, n, &alpha, x, incx, ap);
#		else
			cblas_zspr(ord, ul, n, &alpha, x, incx, ap);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(ap);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SPR_H */

/* EOF */