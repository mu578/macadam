//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_symm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?symm performs the matrix-matrix operation:
 *    c=alpha*a*b + beta*c or c=alpha*b*a + beta*c.
 *
 * \synopsis
 *    void ?symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 *    real-floating alpha, beta
 *    int           lda, ldb, ldc, m, n
 *    char          side, uplo
 *    real-floating a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?symm performs the matrix-matrix operation: c=alpha*a*b + beta*c or or c=alpha*b*a + beta*c
 *    where alpha and beta are scalars, `a` is a symmetric matrix and `b` and `c` are m by n matrices.
 *
 * \parameters
 *    [in] side  - char. Specifies whether  the symmetric matrix `a` appears on the left or right
 *    in the  operation as follows:
 *    side='L' or 'l' c=alpha*a*b + beta*c.
 *    side='R' or 'r' c=alpha*a*b + beta*c.
 *
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the symmetric
 *    matrix `a` is to be referenced as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] m     - int. Specifies the number of rows of the matrix `c`, m must be at least zero.
 *    [in] n     - int. Specifies the number of columns of the matrix `c`, n must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, ka), where ka is m when side='L' or 'l' and is n otherwise.
 *    The m by m part of the array `a` must contain the symmetric matrix, such that when uplo='U' or 'u', the leading
 *    m by m upper triangular part of the array `a` must contain the upper triangular part of the symmetric matrix and
 *    the strictly lower triangular part of `a` is not referenced, and when uplo='L' or 'l', the leading m by m lower
 *    triangular part of the array `a` must contain the lower triangular part of the symmetric matrix and the strictly
 *    upper triangular part of `a` is not referenced.
 *
 *    With side='R' or 'r', the n by n part of the array `a` must contain the symmetric matrix, such that when
 *    uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced, and when
 *    uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When side='L' or 'l' then lda must be at least max(1, m),
 *    otherwise lda must be at least max(1, n).
 *
 *    [int] b    - real-floating array of dimension (ldb, n).
 *    The leading m by n part of the array `b` must contain the matrix `b`.
 *
 *    [in] ldb   - int. Specifies the first dimension of `b`. ldb must be at least max(1, m).
 *
 *    [in] beta  - real-floating. Specifies the scalar beta. When beta is zero then `c` need not be set on input.

 *    [out] c    - real-floating array of dimension (ldc, n).
 *    The leading m by n part of the array `c` must contain the matrix `c`, except when beta is zero. The array `c` is
 *    overwritten by the m by n updated resulted matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
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

#ifndef MC_BLAS_NATIVE_SYMM_H
#define MC_BLAS_NATIVE_SYMM_H

#pragma mark - mc_blas_native_ssymm -

MC_TARGET_FUNC void mc_blas_native_ssymm(const char side, const char uplo, const int m, const int n, const float alpha, const float * a, const int lda, const float * b, const int ldb, const float beta, float * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_SIDE sa = mc_blas_lsame(side, 'L') ? CblasLeft  : CblasRight;
	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_ssymm(ord, sa, ul, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	cblas_ssymm(ord, sa, ul, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#	endif
}

#pragma mark - mc_blas_native_dsymm -

MC_TARGET_FUNC void mc_blas_native_dsymm(const char side, const char uplo, const int m, const int n, const double alpha, const double * a, const int lda, const double * b, const int ldb, const double beta, double * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_SIDE sa = mc_blas_lsame(side, 'L') ? CblasLeft  : CblasRight;
	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_dsymm(ord, sa, ul, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	cblas_dsymm(ord, sa, ul, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#	endif
}

/* \name
 *    ?symm performs the matrix-matrix operation:
 *    c=alpha*a*b + beta*c or c=alpha*b*a + beta*c.
 *
 * \synopsis
 *    void ?symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 *    complex alpha, beta
 *    int     lda, ldb, ldc, m, n
 *    char    side, uplo
 *    complex a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?symm performs the matrix-matrix operation: c=alpha*a*b + beta*c or or c=alpha*b*a + beta*c
 *    where alpha and beta are scalars, `a` is a symmetric matrix and `b` and `c` are m by n matrices.
 *
 * \parameters
 *    [in] side  - char. Specifies whether  the symmetric matrix `a` appears on the left or right
 *    in the  operation as follows:
 *    side='L' or 'l' c=alpha*a*b + beta*c.
 *    side='R' or 'r' c=alpha*a*b + beta*c.
 *
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the symmetric
 *    matrix `a` is to be referenced as follows:
 *    uplo='U' or 'u', the upper triangular part of `a` is being supplied.
 *    uplo='L' or 'l', the lower triangular part of `a` is being supplied.
 *
 *    [in] m     - int. Specifies the number of rows of the matrix `c`, m must be at least zero.
 *    [in] n     - int. Specifies the number of columns of the matrix `c`, n must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, ka), where ka is m when side='L' or 'l' and is n otherwise.
 *    The m by m part of the array `a` must contain the symmetric matrix, such that when uplo='U' or 'u', the leading
 *    m by m upper triangular part of the array `a` must contain the upper triangular part of the symmetric matrix and
 *    the strictly lower triangular part of `a` is not referenced, and when uplo='L' or 'l', the leading m by m lower
 *    triangular part of the array `a` must contain the lower triangular part of the symmetric matrix and the strictly
 *    upper triangular part of `a` is not referenced.
 *
 *    With side='R' or 'r', the n by n part of the array `a` must contain the symmetric matrix, such that when
 *    uplo='U' or 'u', the leading n by n upper triangular part of the array `a` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `a` is not referenced, and when
 *    uplo='L' or 'l', the leading n by n lower triangular part of the array `a` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `a` is not referenced.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When side='L' or 'l' then lda must be at least max(1, m),
 *    otherwise lda must be at least max(1, n).
 *
 *    [int] b    - complex array of dimension (lda, n).
 *    The leading m by n part of the array `b` must contain the matrix `b`.
 *
 *    [in] ldb   - int. Specifies the first dimension of `b`. ldb must be at least max(1, m).
 *
 *    [in] beta  - complex. Specifies the scalar beta. When beta is zero then `c` need not be set on input.

 *    [out] c    - complex array of dimension (ldc, n).
 *    The leading m by n part of the array `c` must contain the matrix `c`, except when beta is zero. The array `c` is
 *    overwritten by the m by n updated resulted matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Jeremy Du Croz, Nag Central Office.
 *     \author Sven Hammarling, Nag Central Office.
 *     \author Richard Hanson, Sandia National Labs.
 */

#pragma mark - mc_blas_native_csymm -

MC_TARGET_FUNC void mc_blas_native_csymm(const char side, const char uplo, const int m, const int n, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * b, const int ldb, const mc_complex_float_t beta, mc_complex_float_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_SIDE sa = mc_blas_lsame(side, 'L') ? CblasLeft  : CblasRight;
	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_csymm(ord, sa, ul, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	else
	cblas_csymm(ord, sa, ul, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	endif
}

#pragma mark - mc_blas_native_zsymm -

MC_TARGET_FUNC void mc_blas_native_zsymm(const char side, const char uplo, const int m, const int n, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * b, const int ldb, const mc_complex_double_t beta, mc_complex_double_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_SIDE sa = mc_blas_lsame(side, 'L') ? CblasLeft  : CblasRight;
	const enum CBLAS_UPLO ul = mc_blas_lsame(uplo, 'U') ? CblasUpper : CblasLower;

#	if MC_TARGET_CPP98
	::cblas_zsymm(ord, sa, ul, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	else
	cblas_zsymm(ord, sa, ul, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SYMM_H */

/* EOF */