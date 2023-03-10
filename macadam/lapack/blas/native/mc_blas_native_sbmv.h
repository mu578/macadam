//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_sbmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?sbmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?sbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
 *    real-floating alpha, beta
 *    int           incx, incy, k, lda, n
 *    char          uplo
 *    real-floating a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?sbmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `a` is an n by n symmetric band matrix, with k
 *    super-diagonals. It computes the matrix-vector product for a real symmetric band matrix,
 *    The band matrix A is stored in either upper or lower-band-packed storage mode, it uses
 *    the scalars alpha and beta, vectors x and y, and band matrix `a`.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the band matrix `a`
 *    is being supplied as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] n     - int. Specifies the ord of the symmetric matrix `a`, n must be at least zero.
 *    [in] k     - int. Specifies the number of super-diagonals of the matrix symmetric `a`, k
 *    must satisfy 0 < k, i.e must be at least one.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row (k + 1) of the array, the first super-diagonal starting at position 1 in row k, and so on. The
 *    top left k by k triangle of the array `a` is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the lower triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row 0 of the array, the first sub-diagonal starting at position 0 in row 1, and so on. The bottom
 *    right k by k triangle of the array `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. `a` must be at least (k + 1).
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
 *              | 1 | 1 | 1 | 1 | 0 | 0 | 0 |
 *              | 1 | 2 | 2 | 2 | 2 | 0 | 0 |
 *              | 1 | 2 | 3 | 3 | 3 | 3 | 0 |
 *     a[7x7] = | 1 | 2 | 3 | 4 | 4 | 4 | 4 |
 *              | 0 | 2 | 3 | 4 | 5 | 5 | 5 |
 *              | 0 | 0 | 3 | 4 | 5 | 6 | 6 |
 *              | 0 | 0 | 0 | 4 | 5 | 6 | 7 |
 *
 *     const real-floating a_band[] = {
 *          0, 0, 0, 1, 2, 3, 4
 *        , 0, 0, 1, 2, 3, 4, 5
 *        , 0, 1, 2, 3, 4, 5, 6
 *        , 1, 2, 3, 4, 5, 6, 7
 *        , 0, 0, 0, 0, 0, 0, 0
 *     };
 *     const real-floating x[] = { 1, 2, 3, 4, 5, 6, 7 };
 *           real-floating y[] = { 1, 0 , 2, 0 , 3, 0 , 4, 0 , 5, 0 , 6, 0 , 7 };
 *     mc_blas_?sbmv('U', 7, 3, 2, a_band, 5, x, 1, 10, y, 2);
 *     on output -> y = { 30, 0 , 78, 0 , 148, 0 , 244, 0 , 288, 0 , 316, 0 , 322 }
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

#ifndef MC_BLAS_NATIVE_SBMV_H
#define MC_BLAS_NATIVE_SBMV_H

#pragma mark - mc_blas_native_ssbmv -

MC_TARGET_FUNC void mc_blas_native_ssbmv(const char uplo, const int n, const int k, const float alpha, const float * a, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_ssbmv(ord, ul, n, k, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_ssbmv(ord, ul, n, k, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

#pragma mark - mc_blas_native_dsbmv -

MC_TARGET_FUNC void mc_blas_native_dsbmv(const char uplo, const int n, const int k, const double alpha, const double * a, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_dsbmv(ord, ul, n, k, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_dsbmv(ord, ul, n, k, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

/* \name
 *    ?sbmv performs the matrix-vector operation:
 *    y=alpha*a*x + beta*y
 *
 * \synopsis
 *    void ?sbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
 *    complex alpha, beta
 *    int     incx, incy, k, lda, n
 *    char    uplo
 *    complex a(lda, *), x(*), y(*)
 *
 * \purpose
 *    ?sbmv performs the matrix-vector operation: y=alpha*a*x + beta*y where alpha and beta are
 *    scalars, `x` and `y` are n element vectors and `a` is an n by n symmetric band matrix, with k
 *    super-diagonals. It computes the matrix-vector product for a real symmetric band matrix,
 *    The band matrix A is stored in either upper or lower-band-packed storage mode, it uses
 *    the scalars alpha and beta, vectors x and y, and band matrix `a`.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the band matrix `a`
 *    is being supplied as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] n     - int. Specifies the ord of the symmetric matrix `a`, n must be at least zero.
 *    [in] k     - int. Specifies the number of super-diagonals of the matrix symmetric `a`, k
 *    must satisfy 0 < k, i.e must be at least one.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row (k + 1) of the array, the first super-diagonal starting at position 1 in row k, and so on. The
 *    top left k by k triangle of the array `a` is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the lower triangular
 *    band part of the symmetric matrix, supplied column by column, with the leading diagonal of the matrix
 *    in row 0 of the array, the first sub-diagonal starting at position 0 in row 1, and so on. The bottom
 *    right k by k triangle of the array `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. `a` must be at least (k + 1).
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

#pragma mark - mc_blas_native_csbmv -

MC_TARGET_FUNC void mc_blas_native_csbmv(const char uplo, const int n, const int k, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * x, const int incx, const mc_complex_float_t beta, mc_complex_float_t * y, const int incy)
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
			::cblas_csbmv(ord, ul, n, k, &alpha, a, lda, x, incx, &beta, y, incy);
#		else
			cblas_csbmv(ord, ul, n, k, &alpha, a, lda, x, incx, &beta, y, incy);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(k);
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

#pragma mark - mc_blas_native_zsbmv -

MC_TARGET_FUNC void mc_blas_native_zsbmv(const char uplo, const int n, const int k, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * x, const int incx, const mc_complex_double_t beta, mc_complex_double_t * y, const int incy)
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
			::cblas_zsbmv(ord, ul, n, k, &alpha, a, lda, x, incx, &beta, y, incy);
#		else
			cblas_zsbmv(ord, ul, n, k, &alpha, a, lda, x, incx, &beta, y, incy);
#		endif
#	else
	mc_unused(uplo);
	mc_unused(n);
	mc_unused(k);
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

#endif /* !MC_BLAS_NATIVE_SBMV_H */

/* EOF */