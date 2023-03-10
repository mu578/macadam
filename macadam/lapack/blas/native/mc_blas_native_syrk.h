//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_syrk.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?syrk performs a rank k operation:
 *    c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c.
 *
 * \synopsis
 *    void ?syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
 *    real-floating alpha, beta
 *    int           k, lda, ldc, n
 *    char          trans, uplo
 *    real-floating a(lda, *), c(ldc, *)
 *
 * \purpose
 *    ?syrk performs a rank k operation: c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c
 *    where alpha and beta are scalars, `c` is an n by n symmetric matrix and `a` is n by k matrix
 *    in the first case and a k by n matrix in the second case.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `c`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `c` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `c` is to be referenced.
 * 
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' c=alpha*a*a' + beta*c + beta*c.
 *    trans='T' or 't' c=alpha*a'*a + beta*c.
 *    trans='C' or 'c' c=alpha*a'*a + beta*c.
 *
 *    [in] n     - int. Specifies the ord of the matrix `c`, n must be at least zero.
 *    [in] k     - int. With trans='N' or 'n', k specifies the number of columns of the matrix `a`,
 *    and with trans='T' or 't' or trans='C' or 'c', k specifies the number of rows of the matrix `a`.
 *    k must be at least zero.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, ka), where ka is k when trans='N' or 'n', and is n otherwise.
 *    With trans='T' or 't', the leading n by k part of the array `a` must contain the matrix `a`, otherwise the leading
 *    k by n part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When trans='N' or 'n' then lda must be at least max(1, n),
 *    otherwise lda must be at least max(1, k).
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
#include <macadam/lapack/blas/mc_blas_xerbla.h>
#include <macadam/details/math/mc_cadd.h>
#include <macadam/details/math/mc_ciseq.h>
#include <macadam/details/math/mc_cmul.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_native_SYRK_H
#define MC_BLAS_native_SYRK_H

#pragma mark - mc_blas_native_ssyrk -

MC_TARGET_FUNC void mc_blas_native_ssyrk(const char uplo, const char trans, const int n, const int k, const float alpha, const float * a, const int lda, const float beta, float * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_ssyrk(ord, ul, tc, n, k, alpha, a, lda, beta, c, ldc);
#	else
	cblas_ssyrk(ord, ul, tc, n, k, alpha, a, lda, beta, c, ldc);
#	endif
}

#pragma mark - mc_blas_native_dsyrk -

MC_TARGET_FUNC void mc_blas_native_dsyrk(const char uplo, const char trans, const int n, const int k, const double alpha, const double * a, const int lda, const double beta, double * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_dsyrk(ord, ul, tc, n, k, alpha, a, lda, beta, c, ldc);
#	else
	cblas_dsyrk(ord, ul, tc, n, k, alpha, a, lda, beta, c, ldc);
#	endif
}

/* \name
 *    ?syrk performs a rank k operation:
 *    c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c.
 *
 * \synopsis
 *    void ?syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
 *    complex alpha, beta
 *    int     k, lda, ldc, n
 *    char    trans, uplo
 *    complex a(lda, *), c(ldc, *)
 *
 * \purpose
 *    ?syrk performs a rank k operation: c=alpha*a*a' + beta*c + beta*c or c=alpha*a'*a + beta*c
 *    where alpha and beta are scalars, `c` is an n by n symmetric matrix and `a` is n by k matrix
 *    in the first case and a k by n matrix in the second case.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the upper or lower triangular part of the array `c`
 *    is to be referenced as follows:
 *    uplo='U' or 'u', only the upper triangular part of `c` is to be referenced.
 *    uplo='L' or 'l', only the lower triangular part of `c` is to be referenced.
 * 
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' c=alpha*a*a' + beta*c + beta*c.
 *    trans='T' or 't' c=alpha*a'*a + beta*c.
 *
 *    [in] n     - int. Specifies the ord of the matrix `c`, n must be at least zero.
 *    [in] k     - int. With trans='N' or 'n', k specifies the number of columns of the matrix `a`,
 *    and with trans='T' or 't' or trans='C' or 'c', k specifies the number of rows of the matrix `a`.
 *    k must be at least zero.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, ka), where ka is k when trans='N' or 'n', and is n otherwise.
 *    With trans='T' or 't', the leading n by k part of the array `a` must contain the matrix `a`, otherwise the
 *    leading k by n part of the array `a` must contain the matrix `a`.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. When trans='N' or 'n' then lda must be at least max(1, n),
 *    otherwise lda must be at least max(1, k).
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

#pragma mark - mc_blas_native_csyrk -

MC_TARGET_FUNC void mc_blas_native_csyrk(const char uplo, const char trans, const int n, const int k, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t beta, mc_complex_float_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_csyrk(ord, ul, tc, n, k, &alpha, a, lda, &beta, c, ldc);
#	else
	cblas_csyrk(ord, ul, tc, n, k, &alpha, a, lda, &beta, c, ldc);
#	endif
}

#pragma mark - mc_blas_native_zsyrk -

MC_TARGET_FUNC void mc_blas_native_zsyrk(const char uplo, const char trans, const int n, const int k, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t beta, mc_complex_double_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE tc = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_zsyrk(ord, ul, tc, n, k, &alpha, a, lda, &beta, c, ldc);
#	else
	cblas_zsyrk(ord, ul, tc, n, k, &alpha, a, lda, &beta, c, ldc);
#	endif
}

#endif /* !MC_BLAS_SYRK_H */

/* EOF */