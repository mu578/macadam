//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_iamax.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    i?amax returns the index of the first element having maximum absolute value.
 *
 * \synopsis
 *    int i?amax(n, x, incx)
 *    int           incx, n
 *    real-floating x(*)
 *
 * \purpose
 *    ?nrm2 returns the index of the first element having maximum absolute value.
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vector `x`.
 *    [in] x     - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Sven Hammarling, Nag Ltd.
 */

#include <macadam/details/mc_target.h>

#ifndef MC_BLAS_IAMAX_H
#define MC_BLAS_IAMAX_H

#pragma mark - mc_blas_native_isamax -

MC_TARGET_FUNC int mc_blas_native_isamax(const int n, const float * x, const int incx)
{
	int iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
#	if MC_TARGET_CPP98
	iamax = ::cblas_isamax(n, x, incx);
#	else
	iamax = cblas_isamax(n, x, incx);
#	endif
	return iamax;
}

#pragma mark - mc_blas_native_idamax -

MC_TARGET_FUNC int mc_blas_native_idamax(const int n, const double * x, const int incx)
{
	int iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
#	if MC_TARGET_CPP98
	iamax = ::cblas_idamax(n, x, incx);
#	else
	iamax = cblas_idamax(n, x, incx);
#	endif
	return iamax;
}

/* \name
 *    i?amax returns the index of the first element having maximum absolute value.
 *
 * \synopsis
 *    int i?amax(n, x, incx)
 *    int     incx, n
 *    complex x(*)
 *
 * \purpose
 *    ?nrm2 returns the index of the first element having maximum absolute value.
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vector `x`.
 *    [in] x     - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#pragma mark - mc_blas_native_icamax -

MC_TARGET_FUNC int mc_blas_native_icamax(const int n, const mc_complex_float_t * x, const int incx)
{
	int iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
#	if MC_TARGET_CPP98
	iamax = ::cblas_icamax(n, x, incx);
#	else
	iamax = cblas_icamax(n, x, incx);
#	endif
	return iamax;
}

#pragma mark - mc_blas_native_izamax -

MC_TARGET_FUNC int mc_blas_native_izamax(const int n, const mc_complex_double_t * x, const int incx)
{
	int iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
#	if MC_TARGET_CPP98
	iamax = ::cblas_izamax(n, x, incx);
#	else
	iamax = cblas_izamax(n, x, incx);
#	endif
	return iamax;
}

#endif /* !MC_BLAS_IAMAX_H */

/* EOF */