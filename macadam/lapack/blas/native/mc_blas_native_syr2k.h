//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_syr2k.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?syr2k performs a rank 2k operation:
 *    c=alpha*a*b' + alpha*b*a' + beta*c or c=alpha*a'*b + alpha*b'*a + beta*c.
 *
 * \synopsis
 *    void ?syr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 *    real-floating alpha, beta
 *    int           k, lda, ldb, ldc, n
 *    char          trans, uplo
 *    real-floating a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?syr2k performs a rank 2k operation: c=alpha*a*b' + alpha*b*a' + beta*c or c=alpha*a'*b + alpha*b'*a + beta*c
 *    where alpha and beta are scalars, `c` is an n by n symmetric matrix and `a` and `b` are n by k matrices in the
 *    first case and k by n matrices in the second case.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `c`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `c` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `c` is to be referenced.
 * 
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' c=alpha*a*b' + alpha*b*a' + beta*c.
 *    trans='T' or 't' c=alpha*a'*b + alpha*b'*a + beta*c.
 *    trans='C' or 'c' c=alpha*a'*b + alpha*b'*a + beta*c.
 *
 *    [in] n     - int. Specifies the ord of the matrix `c`, n must be at least zero.
 *    [in] k     - int. With trans='N' or 'n', k specifies the number of columns of the matrices `a` and `b`,
 *    and with trans='T' or 't' or trans='C' or 'c', k specifies the number of rows of the matrices `a` and `b`.
 *    k must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, ka), where ka is k when trans='N' or 'n', and is n otherwise.
 *    With trans='T' or 't', the leading n by k part of the array `a` must contain the matrix `a`, otherwise the leading
 *    k by n  part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When trans='N' or 'n' then lda must be at least max(1, n),
 *    otherwise lda must be at least max(1, k).
 *
 *    [int] b    - real-floating array of dimension (ldb, kb), where kb is k when trans='N' or 'n', and is n otherwise.
 *    With trans='N' or 'n', the leading n by k part of the array `b` must contain the matrix `b`, otherwise the leading
 *    k by n part of the array `b` must contain the matrix `b`.
 *
 *    [in] ldb   - int. Specifies the first dimension of `b`. When trans='N' or 'n'
 *    then ldb must be at least max(1, n), otherwise ldb must be at least max(1, k).
 *
 *    [in] beta  - real-floating. Specifies the scalar beta. When beta is zero then `c` need not be set on input.
 *
 *    [out] c    - real-floating array of dimension (ldc, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `c` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `c` is not referenced. The upper triangular
 *    part of the array `c` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `c` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `c` is not referenced. The lower triangular
 *    part of the array `c` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, n).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Iain S. Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>

#ifndef MC_BLAS_NATIVE_SYR2K_H
#define MC_BLAS_NATIVE_SYR2K_H

#pragma mark - mc_blas_native_ssyr2k -

MC_TARGET_FUNC void mc_blas_native_ssyr2k(const char uplo, const char trans, const int n, const int k, const float alpha, const float * a, const int lda, const float * b, const int ldb, const float beta, float * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_ssyr2k(ord, ul, tc, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	cblas_ssyr2k(ord, ul, tc, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#	endif
}

#pragma mark - mc_blas_native_dsyr2k -

MC_TARGET_FUNC void mc_blas_native_dsyr2k(const char uplo, const char trans, const int n, const int k, const double alpha, const double * a, const int lda, const double * b, const int ldb, const double beta, double * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_dsyr2k(ord, ul, tc, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#	else
	cblas_dsyr2k(ord, ul, tc, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#	endif
}

/* \name
 *    ?syr2k performs a rank 2k operation:
 *    c=alpha*a*b' + alpha*b*a' + beta*c or c=alpha*a'*b + alpha*b'*a + beta*c.
 *
 * \synopsis
 *    void ?syr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 *    complex alpha, beta
 *    int     k, lda, ldb, ldc, n
 *    char    trans, uplo
 *    complex a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?syr2k performs a rank 2k operation: c=alpha*a*b' + alpha*b*a' + beta*c or c=alpha*a'*b + alpha*b'*a + beta*c
 *    where alpha and beta are scalars, `c` is an n by n symmetric matrix and `a` and `b` are n by k matrices in the
 *    first case and k by n matrices in the second case.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `c`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `c` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `c` is to be referenced.
 * 
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' c=alpha*a*b' + alpha*b*a' + beta*c.
 *    trans='T' or 't' c=alpha*a'*b + alpha*b'*a + beta*c.
 *
 *    [in] n     - int. Specifies the ord of the matrix `c`, n must be at least zero.
 *    [in] k     - int. With trans='N' or 'n', k specifies the number of columns of the matrices `a` and `b`,
 *    and with trans='T' or 't' or trans='C' or 'c', k specifies the number of rows of the matrices `a` and `b`.
 *    k must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, ka), where ka is k when trans='N' or 'n', and is n otherwise.
 *    With trans='T' or 't', the leading n by k part of the array `a` must contain the matrix `a`, otherwise the
 *    leading k by n  part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When trans='N' or 'n' then lda must be at least max(1, n),
 *    otherwise lda must be at least max(1, k).
 *
 *    [int] b    - complex array of dimension (ldb, kb), where kb is k when trans='N' or 'n', and is n otherwise.
 *    With trans='N' or 'n', the leading n by k part of the array `b` must contain the matrix `b`, otherwise the
 *    leading k by n part of the array `b` must contain the matrix `b`.
 *
 *    [in] ldb   - int. Specifies the first dimension of `b`. When trans='N' or 'n'
 *    then ldb must be at least max(1, n), otherwise ldb must be at least max(1, k).
 *
 *    [in] beta  - complex. Specifies the scalar beta. When beta is zero then `c` need not be set on input.
 *
 *    [out] c    - complex array of dimension (ldc, n).
 *    With uplo='U' or 'u', the leading n by n upper triangular part of the array `c` must contain the upper triangular
 *    part of the symmetric matrix and the strictly lower triangular part of `c` is not referenced. The upper triangular
 *    part of the array `c` is overwritten by the upper triangular part of the updated matrix.
 *
 *    With uplo='L' or 'l', the leading n by n lower triangular part of the array `c` must contain the lower triangular
 *    part of the symmetric matrix and the strictly upper triangular part of `c` is not referenced. The lower triangular
 *    part of the array `c` is overwritten by the lower triangular part of the updated matrix.
 *
 *    [in] ldc   - int. Specifies the first dimension of `c`. ldc must be at least max(1, n).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Iain S. Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#pragma mark - mc_blas_native_csyr2k -

MC_TARGET_FUNC void mc_blas_native_csyr2k(const char uplo, const char trans, const int n, const int k, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * b, const int ldb, const mc_complex_float_t beta, mc_complex_float_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_csyr2k(ord, ul, tc, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	else
	cblas_csyr2k(ord, ul, tc, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	endif
}

#pragma mark - mc_blas_native_zsyr2k -

MC_TARGET_FUNC void mc_blas_native_zsyr2k(const char uplo, const char trans, const int n, const int k, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * b, const int ldb, const mc_complex_double_t beta, mc_complex_double_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_zsyr2k(ord, ul, tc, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	else
	cblas_zsyr2k(ord, ul, tc, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SYR2K_H */

/* EOF */