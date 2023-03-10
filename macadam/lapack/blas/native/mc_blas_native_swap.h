//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_swap.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?swap interchanges two vectors.
 *
 * \synopsis
 *    void ?swap(n, x, incx, y, incy)
 *    int           incx, incy, n
 *    real-floating x(*), y(*)
 *
 * \purpose
 *    ?swap interchanges two vectors.
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in input vector(s).
 *
 *    [out] x    - real-floating array of dimension (at least) (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the storage spacing between elements of `x`.
 *
 *    [out]  y   - real-floating array of dimension (at least) (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the storage spacing between elements of `y`.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#include <macadam/lapack/blas/mc_blas_access.h>

#ifndef MC_BLAS_NATIVE_SWAP_H
#define MC_BLAS_NATIVE_SWAP_H

#pragma mark - mc_blas_native_sswap -

MC_TARGET_FUNC void mc_blas_native_sswap(const int n, float * x, const int incx, float * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_sswap(n, x, incx, y, incy);
#	else
	cblas_sswap(n, x, incx, y, incy);
#	endif
}

#pragma mark - mc_blas_native_dswap -

MC_TARGET_FUNC void mc_blas_native_dswap(const int n, double * x, const int incx, double * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_dswap(n, x, incx, y, incy);
#	else
	cblas_dswap(n, x, incx, y, incy);
#	endif
}

/* \name
 *    ?swap interchanges two vectors.
 *
 * \synopsis
 *    void ?swap(n, x, incx, y, incy)
 *    int     incx, incy, n
 *    complex x(*), y(*)
 *
 * \purpose
 *    ?swap interchanges two vectors.
 *
 * \parameters
 *    [in]  n    - int. Specifies the number of elements in input vector(s).
 *
 *    [out] x    - complex array of dimension (at least) (1+(n-1)*abs(incx)).
 *    [in]  incx - int. Specifies the storage spacing between elements of `x`.
 *
 *    [out]  y   - complex array of dimension (at least) (1+(n-1)*abs(incy)).
 *    [in]  incy - int. Specifies the storage spacing between elements of `y`.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#pragma mark - mc_blas_native_cswap -

MC_TARGET_FUNC void mc_blas_native_cswap(const int n, mc_complex_float_t * x, const int incx, mc_complex_float_t * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_cswap(n, x, incx, y, incy);
#	else
	cblas_cswap(n, x, incx, y, incy);
#	endif
}

#pragma mark - mc_blas_native_zswap -

MC_TARGET_FUNC void mc_blas_native_zswap(const int n, mc_complex_double_t * x, const int incx, mc_complex_double_t * y, const int incy)
{
#	if MC_TARGET_CPP98
	::cblas_zswap(n, x, incx, y, incy);
#	else
	cblas_zswap(n, x, incx, y, incy);
#	endif
}

#endif /* !MC_BLAS_NATIVE_SWAP_H */

/* EOF */