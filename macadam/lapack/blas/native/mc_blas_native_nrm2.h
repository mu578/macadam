//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_nrm2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?nrm2 returns the euclidean norm of a vector.
 *
 * \synopsis
 *    real-floating ?nrm2(n, x, incx)
 *    int           incx, n
 *    real-floating x(*)
 *
 * \purpose
 *    ?nrm2 returns the euclidean norm of a vector: norm2=sqrt(x'*x).
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

#include <macadam/lapack/blas/mc_blas_access.h>

#ifndef MC_BLAS_NATIVE_NRM2_H
#define MC_BLAS_NATIVE_NRM2_H

#pragma mark - mc_blas_native_snrm2 -

MC_TARGET_FUNC float mc_blas_native_snrm2(const int n, const float * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_snrm2(n, x, incx);
#	else
	return cblas_snrm2(n, x, incx);
#	endif
}

#pragma mark - mc_blas_native_dnrm2 -

MC_TARGET_FUNC double mc_blas_native_dnrm2(const int n, const double * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_dnrm2(n, x, incx);
#	else
	return cblas_dnrm2(n, x, incx);
#	endif
}

/* \name
 *    ?nrm2 returns the euclidean norm of a vector.
 *
 * \synopsis
 *    real-floating ?nrm2(n, x, incx)
 *    int     incx, n
 *    complex x(*)
 *
 * \purpose
 *    ?nrm2 returns the euclidean norm of a vector: norm2=sqrt(x_*x).
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
 *     \author Sven Hammarling, Nag Ltd.
 */

#pragma mark - mc_blas_native_scnrm2 -

MC_TARGET_FUNC float mc_blas_native_scnrm2(const int n, const mc_complex_float_t * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_scnrm2(n, x, incx);
#	else
	return cblas_scnrm2(n, x, incx);
#	endif
}

#pragma mark - mc_blas_native_dznrm2 -

MC_TARGET_FUNC double mc_blas_native_dznrm2(const int n, const mc_complex_double_t * x, const int incx)
{
#	if MC_TARGET_CPP98
	return ::cblas_dznrm2(n, x, incx);
#	else
	return cblas_dznrm2(n, x, incx);
#	endif
}

#endif /* !MC_BLAS_NATIVE_NRM2_H */

/* EOF */