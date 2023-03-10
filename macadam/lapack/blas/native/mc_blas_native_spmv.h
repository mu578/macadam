//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_spmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?spmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?spmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
 *    real-floating alpha, beta
 *    int           incx, incy, n
 *    char          uplo
 *    real-floating ap(*), x(*), y(*)
 *
 * \purpose
 *    ?spmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `ap` is an n by n symmetric matrix, supplied in
 *    packed form of dimension (at least) ((n*(n+1))/2).
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the matrix `a` is
 *    supplied in the packed array `ap` as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` supplied in `ap`.
 *    uplo='L' or 'l', the lower triangular part of `a` supplied in `ap`.
 *
 *    [in] n     - int. Specifies the ord of the symmetric matrix `a`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
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
 *    [int] x    - real-floating array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must
 *    contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - real-floating. Specifies the scalar beta.
 *
 *    [out] y    - real-floating array of size at least (1+(n-1)*abs(incy)). The incremented array `y` must
 *    contain the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
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

#ifndef MC_BLAS_NATIVE_SPMV_H
#define MC_BLAS_NATIVE_SPMV_H

#pragma mark - mc_blas_native_sspmv -

MC_TARGET_FUNC void mc_blas_native_sspmv(const char uplo, const int n, const float alpha, const float * ap, const float * x, const int incx, const float beta, float * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_sspmv(ord, ul, n, alpha, ap, x, incx, beta, y, incy);
#	else
	cblas_sspmv(ord, ul, n, alpha, ap, x, incx, beta, y, incy);
#	endif
}

#pragma mark - mc_blas_native_dspmv -

MC_TARGET_FUNC void mc_blas_native_dspmv(const char uplo, const int n, const double alpha, const double * ap, const double * x, const int incx, const double beta, double * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_dspmv(ord, ul, n, alpha, ap, x, incx, beta, y, incy);
#	else
	cblas_dspmv(ord, ul, n, alpha, ap, x, incx, beta, y, incy);
#	endif
}

/* \name
 *    ?spmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?spmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
 *    complex alpha, beta
 *    int     incx, incy, n
 *    char    uplo
 *    complex ap(*), x(*), y(*)
 *
 * \purpose
 *    ?spmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `ap` is an n by n symmetric matrix, supplied in
 *    packed form of dimension (at least) ((n*(n+1))/2).
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the matrix `a` is
 *    supplied in the packed array `ap` as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` supplied in `ap`.
 *    uplo='L' or 'l', the lower triangular part of `a` supplied in `ap`.
 *
 *    [in] n     - int. Specifies the ord of the symmetric matrix `a`, n must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] ap    - complex array of dimension (at least) ((n*(n+1))/2).
 *    With uplo='U' or 'u', the array `ap` must contain the upper triangular part of the symmetric matrix
 *    packed sequentially, if column-major layout: column by column, so that ap(1) contains a(1,1), ap(2)
 *    and ap(3) contain a(1,2) and a(2,2) respectively, and so on, else row by row following the same logic.
 *
 *    With uplo='L' or 'l', the array `a`p must contain the lower triangular part of the symmetric matrix
 *    packed sequentially, if column-major layout: column by column, so that ap(1) contains a(1,1), ap(2)
 *    and ap(3) contain a(2,1)  and a(3,1) respectively, and so on, else row by row following the same logic.
 *
 *    [int] x    - complex array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must
 *    contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - complex. Specifies the scalar beta.
 *
 *    [out] y    - complex array of size at least (1+(n-1)*abs(incy)). The incremented array `y` must
 *    contain the vector `y`, y is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
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

#pragma mark - mc_blas_native_cspmv -

MC_TARGET_FUNC void mc_blas_native_cspmv(const char uplo, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * ap, const mc_complex_float_t * x, const int incx, const mc_complex_float_t beta, mc_complex_float_t * y, const int incy)
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
			::cblas_cspmv(ord, ul, n, &alpha, ap, x, incx, &beta, y, incy);
#		else
			cblas_cspmv(ord, ul, n, &alpha, ap, x, incx, &beta, y, incy);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(ap);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(beta);
	mc_unused(y);
	mc_unused(incy);
#	endif
}

#pragma mark - mc_blas_native_zspmv -

MC_TARGET_FUNC void mc_blas_native_zspmv(const char uplo, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * ap, const mc_complex_double_t * x, const int incx, const mc_complex_double_t beta, mc_complex_double_t * y, const int incy)
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
			::cblas_zspmv(ord, ul, n, &alpha, ap, x, incx, &beta, y, incy);
#		else
			cblas_zspmv(ord, ul, n, &alpha, ap, x, incx, &beta, y, incy);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(alpha);
	mc_unused(ap);
	mc_unused(x);
	mc_unused(incx);
	mc_unused(beta);
	mc_unused(y);
	mc_unused(incy);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SPMV_H */

/* EOF */