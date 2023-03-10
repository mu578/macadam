//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_gemm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?gemm performs one of the matrix-matrix operations:
 *    c=alpha*op(a)*op(b) + beta*c where op(x)=x or op(x)=x'.
 *
 * \synopsis
 *    void ?gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 *    real-floating alpha, beta
 *    int           k, lda, ldb, ldc, m, n
 *    char          transa, transb
 *    real-floating a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?gemm performs one of the matrix-matrix operations: c=alpha*op(a)*op(b) + beta*c where
 *    op(x)=x or op(x)=x' alpha and beta are scalars, and a, b and c are matrices, with op(a)
 *    an m by k matrix, op(b) a k by n matrix and c an m by n matrix.
 *
 * \parameters
 *    [in] transa - char. Specifies the form of op(a) to be used in the matrix multiplication as follows:
 *    transa='N' or 'n', op(a)=a.
 *    transa='T' or 't', op(a)=a'.
 *    transa='C' or 'c', op(a)=a'.
 *
 *    [in] transb - char. Specifies the form of op(b) to be used in the matrix multiplication as follows:
 *    transb='N' or 'n', op(b)=b.
 *    transb='T' or 't', op(b)=b'.
 *    transb='C' or 'c', op(b)=b'.
 *
 *    [in] m      - int. Specifies the number of rows of the matrix op(a) and of the matrix `c`, m must be
 *    at least zero.
 *
 *    [in] n      - int. Specifies the number of columns of the matrix op(b) and the number of columns of
 *    the matrix `c`, n must be at least zero.
 *
 *    [in] k      - int. Specifies  the number of columns of the matrix op(a) and the number of rows of
 *    the matrix op(b), k must be at least zero.
 *
 *    [in] alpha  - real-floating. Specifies the scalar alpha.
 *
 *    [in] a      - real-floating array of dimension (lda, ka), where ka is k when transa='N' or 'n' and
 *    is m otherwise. Prior entry with transa='N' or 'n', the leading m by k part of the array a must
 *    contain the matrix `a`, otherwise the leading k by m part of the array a must contain the matrix `a`.
 *
 *    [in] lda    - int. Specifies the first dimension of `a`. When transa='N' or 'n' then
 *    lda must be at least max(1, m), otherwise lda must be at least max(1, k).
 *
 *    [in] b      - real-floating array of dimension (ldb, kb), where kb is n when transb='N' or 'n' and
 *    is k otherwise. Prior entry with transb='N' or 'n', the leading k by n part of the array b must
 *    contain the matrix `b`, otherwise the leading n by k part of the array b must contain the matrix `b`.
 *
 *    [in] ldb    - int. Specifies the first dimension of `b`. When transb='N' or 'n' then
 *    ldb must be at least max(1, k), otherwise ldb must be at least max(1, n).
 *
 *    [in] beta   - real-floating. Specifies the scalar beta. When beta is supplied as zero then c need
 *    not be set on input.
 *
 *    [out] c     - real-floating array of dimension (ldc, n). Prior entry the leading  m by n part of the
 *    array c must contain the matrix `c`, except when beta is set to zero, in which case c need not be set
 *    on entry, c is overwritten by the m by n matrix (alpha*op(a)*op(b) + beta*c).
 *
 *    [in] ldc    - int. Specifies the first dimension of `c`, ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Laboratory.
 *     \author Iain Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#include <macadam/lapack/blas/mc_blas_lsame.h>

#ifndef MC_BLAS_NATIVE_GEMM_H
#define MC_BLAS_NATIVE_GEMM_H

#pragma mark - mc_blas_native_sgemm -

MC_TARGET_FUNC void mc_blas_native_sgemm(const char transa, const char transb, const int m, const int n, const int k, const float alpha, const float * a, const int lda, const float * b, const int ldb, const float beta, float * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
	const int ma = n;
	const int mb = k;
	const int mc = n;
	mc_unused(lda);
	mc_unused(ldb);
	mc_unused(ldc);
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
	const int ma = lda;
	const int mb = ldb;
	const int mc = ldc;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(transa, 'N') ? CblasNoTrans : (mc_blas_lsame(transa, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_TRANSPOSE tb = mc_blas_lsame(transb, 'N') ? CblasNoTrans : (mc_blas_lsame(transb, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_sgemm(ord, ta, tb, m, n, k, alpha, a, ma, b, mb, beta, c, mc);
#	else
	cblas_sgemm(ord, ta, tb, m, n, k, alpha, a, ma, b, mb, beta, c, mc);
#	endif
}

#pragma mark - mc_blas_native_dgemm -

MC_TARGET_FUNC void mc_blas_native_dgemm(const char transa, const char transb, const int m, const int n, const int k, const double alpha, const double * a, const int lda, const double * b, const int ldb, const double beta, double * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
	const int ma = n;
	const int mb = k;
	const int mc = n;
	mc_unused(lda);
	mc_unused(ldb);
	mc_unused(ldc);
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
	const int ma = lda;
	const int mb = ldb;
	const int mc = ldc;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(transa, 'N') ? CblasNoTrans : (mc_blas_lsame(transa, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_TRANSPOSE tb = mc_blas_lsame(transb, 'N') ? CblasNoTrans : (mc_blas_lsame(transb, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_dgemm(ord, ta, tb, m, n, k, alpha, a, ma, b, mb, beta, c, mc);
#	else
	cblas_dgemm(ord, ta, tb, m, n, k, alpha, a, ma, b, mb, beta, c, mc);
#	endif
}

/* \name
 *    ?gemm performs one of the matrix-matrix operations:
 *    c=alpha*op(a)*op(b) + beta*c where op(x)=x or op(x)=x' or op(x)=x_.
 *
 * \synopsis
 *    void ?gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 *    complex alpha, beta
 *    int     k, lda, ldb, ldc, m, n
 *    char    transa, transb
 *    complex a(lda, *), b(ldb, *), c(ldc, *)
 *
 * \purpose
 *    ?gemm performs one of the matrix-matrix operations: c=alpha*op(a)*op(b) + beta*c where
 *    op(x)=x or op(x)=x' or op(x)=x_ alpha and beta are scalars, and a, b and c are matrices,
 *    with op(a) an m by k matrix, op(b) a k by n matrix and c an m by n matrix.
 *
 * \parameters
 *    [in] transa - char. Specifies the form of op(a) to be used in the matrix multiplication as follows:
 *    transa='N' or 'n', op(a)=a.
 *    transa='T' or 't', op(a)=a'.
 *    transa='C' or 'c', op(a)=a_.
 *
 *    [in] transb - char. Specifies the form of op(b) to be used in the matrix multiplication as follows:
 *    transb='N' or 'n', op(b)=b.
 *    transb='T' or 't', op(b)=b'.
 *    transb='C' or 'c', op(b)=b_.
 *
 *    [in] m      - int. Specifies the number of rows of the matrix op(a) and of the matrix `c`, m must be
 *    at least zero.
 *
 *    [in] n      - int. Specifies the number of columns of the matrix op(b) and the number of columns of
 *    the matrix `c`, n must be at least zero.
 *
 *    [in] k      - int. Specifies  the number of columns of the matrix op(a) and the number of rows of
 *    the matrix op(b), k must be at least zero.
 *
 *    [in] alpha  - complex. Specifies the scalar alpha.
 *
 *    [in] a      - complex array of dimension (lda, ka), where ka is k when transa='N' or 'n' and
 *    is m otherwise. Prior entry with transa='N' or 'n', the leading m by k part of the array a must
 *    contain the matrix `a`, otherwise the leading k by m part of the array a must contain the matrix `a`.
 *
 *    [in] lda    - int. Specifies the first dimension of `a`. When transa='N' or 'n' then
 *    lda must be at least max(1, m), otherwise lda must be at least max(1, k).
 *
 *    [in] b      - complex array of dimension (ldb, kb), where kb is n when transb='N' or 'n' and
 *    is k otherwise. Prior entry with transb='N' or 'n', the leading k by n part of the array b must
 *    contain the matrix `b`, otherwise the leading n by k part of the array b must contain the matrix `b`.
 *
 *    [in] ldb    - int. Specifies the first dimension of `b`. When transb='N' or 'n' then
 *    ldb must be at least max(1, k), otherwise ldb must be at least max(1, n).
 *
 *    [in] beta   - complex. Specifies the scalar beta. When beta is supplied as zero then c need
 *    not be set on input.
 *
 *    [out] c     - complex array of dimension (ldc, n). Prior entry the leading  m by n part of the
 *    array c must contain the matrix `c`, except when beta is set to zero, in which case c need not be set
 *    on entry, c is overwritten by the m by n matrix (alpha*op(a)*op(b) + beta*c).
 *
 *    [in] ldc    - int. Specifies the first dimension of `c`, ldc must be at least max(1, m).
 *
 * \examples
 *
 * \level 3 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Laboratory.
 *     \author Iain Duff, AERE Harwell.
 *     \author Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     \author Sven Hammarling, Numerical Algorithms Group Ltd.
 */

#pragma mark - mc_blas_native_cgemm -

MC_TARGET_FUNC void mc_blas_native_cgemm(const char transa, const char transb, const int m, const int n, const int k, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * b, const int ldb, const mc_complex_float_t beta, mc_complex_float_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
	const int ma = n;
	const int mb = k;
	const int mc = n;
	mc_unused(lda);
	mc_unused(ldb);
	mc_unused(ldc);
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
	const int ma = lda;
	const int mb = ldb;
	const int mc = ldc;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(transa, 'N') ? CblasNoTrans : (mc_blas_lsame(transa, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_TRANSPOSE tb = mc_blas_lsame(transb, 'N') ? CblasNoTrans : (mc_blas_lsame(transb, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_cgemm(ord, ta, tb, m, n, k, &alpha, a, ma, b, mb, &beta, c, mc);
#	else
	cblas_cgemm(ord, ta, tb, m, n, k, &alpha, a, ma, b, mb, &beta, c, mc);
#	endif
}

#pragma mark - mc_blas_native_zgemm -

MC_TARGET_FUNC void mc_blas_native_zgemm(const char transa, const char transb, const int m, const int n, const int k, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * b, const int ldb, const mc_complex_double_t beta, mc_complex_double_t * c, const int ldc)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
	const int ma = n;
	const int mb = k;
	const int mc = n;
	mc_unused(lda);
	mc_unused(ldb);
	mc_unused(ldc);
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
	const int ma = lda;
	const int mb = ldb;
	const int mc = ldc;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(transa, 'N') ? CblasNoTrans : (mc_blas_lsame(transa, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_TRANSPOSE tb = mc_blas_lsame(transb, 'N') ? CblasNoTrans : (mc_blas_lsame(transb, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_zgemm(ord, ta, tb, m, n, k, &alpha, a, ma, b, mb, &beta, c, mc);
#	else
	cblas_zgemm(ord, ta, tb, m, n, k, &alpha, a, ma, b, mb, &beta, c, mc);
#	endif
}

#endif /* !MC_BLAS_NATIVE_GEMM_H */

/* EOF */