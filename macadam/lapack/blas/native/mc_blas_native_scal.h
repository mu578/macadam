//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_scal.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?scal scales a vector by a constant.
 *
 * \synopsis
 *    void ?scal(n, a, x, incx)
 *    int           incx, n
 *    real-floating a
 *    real-floating x(*)
 *
 * \purpose
 *    ?scal scales a vector by a constant.
 *
 * \parameters
 *    [in] n    - int. Specifies the number of elements in the input vector `x`.
 *    [in] a    - real-floating. Specifies the scalar alpha.
 *    [in] x    - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
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

#ifndef MC_BLAS_NATIVE_SCAL_H
#define MC_BLAS_NATIVE_SCAL_H

#pragma mark - mc_blas_native_sscal -

MC_TARGET_FUNC void mc_blas_native_sscal(const int n, float a, float * x, const int incx)
{
#	if MC_TARGET_CPP98
	::cblas_sscal(n, a, x, incx);
#	else
	cblas_sscal(n, a, x, incx);
#	endif
}

#pragma mark - mc_blas_native_dscal -

MC_TARGET_FUNC void mc_blas_native_dscal(const int n, double a, double * x, const int incx)
{
#	if MC_TARGET_CPP98
	::cblas_dscal(n, a, x, incx);
#	else
	cblas_dscal(n, a, x, incx);
#	endif
}

/* \name
 *    ?scal scales a vector by a constant.
 *
 * \synopsis
 *    void ?scal(n, a, x, incx)
 *    int     incx, n
 *    complex a
 *    complex x(*)
 *
 * \purpose
 *    ?scal scales a vector by a constant.
 *
 * \parameters
 *    [in] n    - int. Specifies the number of elements in the input vector `x`.
 *    [in] a    - complex. Specifies the scalar alpha.
 *    [in] x    - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#pragma mark - mc_blas_native_cscal -

MC_TARGET_FUNC void mc_blas_native_cscal(const int n, mc_complex_float_t a, mc_complex_float_t * x, const int incx)
{
#	if MC_TARGET_CPP98
	::cblas_cscal(n, &a, x, incx);
#	else
	cblas_cscal(n, &a, x, incx);
#	endif
}

#pragma mark - mc_blas_native_zscal -

MC_TARGET_FUNC void mc_blas_native_zscal(const int n, mc_complex_double_t a, mc_complex_double_t * x, const int incx)
{
#	if MC_TARGET_CPP98
	::cblas_zscal(n, &a, x, incx);
#	else
	cblas_zscal(n, &a, x, incx);
#	endif
}

#endif /* !MC_BLAS_SCAL_H */

/* EOF */