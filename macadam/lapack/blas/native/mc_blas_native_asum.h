//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_asum.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?asum returns the sum of the absolute values of a vector.
 *
 * \synopsis
 *    real-floating ?asum(n, x, incx)
 *    int           incx, n
 *    real-floating x(*)
 *
 * \purpose
 *    ?asum returns the sum of the absolute values of a vector.
 *
 * \parameters
 *    [in] n    - int. Specifies the number of elements in the input vector `x`.
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

#ifndef MC_BLAS_NATIVE_ASUM_H
#define MC_BLAS_NATIVE_ASUM_H

#pragma mark - mc_blas_native_sasum -

MC_TARGET_FUNC float mc_blas_native_sasum(const int n, const float * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_sasum(n, x, incx);
#	else
	return cblas_sasum(n, x, incx);
#	endif
}

#pragma mark - mc_blas_native_dasum -

MC_TARGET_FUNC double mc_blas_native_dasum(const int n, const double * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_dasum(n, x, incx);
#	else
	return cblas_dasum(n, x, incx);
#	endif
}

/* \name
 *    ?asum returns the sum of the absolute values of a vector.
 *
 * \synopsis
 *    real-floating ?asum(n, x, incx)
 *    int     incx, n
 *    complex x(*)
 *
 * \purpose
 *    ?asum returns the sum of the absolute values of a vector.
 *
 * \parameters
 *    [in] n    - int. Specifies the number of elements in the input vector `x`.
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

#pragma mark - mc_blas_native_scasum -

MC_TARGET_FUNC float mc_blas_native_scasum(const int n, const mc_complex_float_t * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_scasum(n, x, incx);
#	else
	return cblas_scasum(n, x, incx);
#	endif
}

#pragma mark - mc_blas_native_dzasum -

MC_TARGET_FUNC double mc_blas_native_dzasum(const int n, const mc_complex_double_t * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_dzasum(n, x, incx);
#	else
	return cblas_dzasum(n, x, incx);
#	endif
}

#endif /* !MC_BLAS_NATIVE_ASUM_H */

/* EOF */