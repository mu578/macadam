//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_rotm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?rotm performs modified Givens rotation of points in the plane.
 *
 * \synopsis
 *    real-floating ?rotm(n, x, incx, y, incy, param)
 *    int           incx, incy, n
 *    real-floating param(5), x(*), y(*)
 *
 * \purpose
 *    ?rotm performs modified Givens rotation of points in the plane. Given two vectors x and y,
 *    each vector element of these vectors is replaced as follows:
 *      | xi |     | xi |
 *      | yi | = H | Yi |
 *
 *    for i=0 to n-1, where H is a modified Givens transformation matrix whose values are stored
 *    in the param[1] through param array.
 *
 *    param(0) contains a switch, flag. param[1 - 4] contain h11, h21, h12, and h22,
 *    respectively, the components of the array H.
 *
 *    Depending on the values of flag, the components of H are set as follows:
 *       flag=-1:
 *       H[2x2] = | h11 | h12 |
 *                | h21 | h22 |
 *       flag=0:
 *       H[2x2] = |  1  | h12 |
 *                | h21 |  1  |
 *       flag=1:
 *       H[2x2] = | h11 |  1  |
 *                | -1  | h22 |
 *       flag=-2:
 *       H[2x2] = |  1  |  0  |
 *                |  0  |  1  |
 *
 *    In the last three cases, the matrix entries of 1, -1, and 0 are assumed based on the
 *    value of flag and are not required to be set in the param array.
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vector `x` and `y`.
 *
 *    [out] x    - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    Each element x[i] is replaced by x[i]=h11*x[i] + h12*y[i].
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [out] y    - real-floating array of size at least (1+(n-1)*abs(incy)).
 *    Each element y[i] is replaced by y[i]=h21*x[i] + h22*y[i].
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`, incy must not be zero.
 *
 *    [in] param - real-floating array of size 5 such as param[flag, h11, h21, h12, h22].
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

#ifndef MC_BLAS_NATIVE_ROTM_H
#define MC_BLAS_NATIVE_ROTM_H

#pragma mark - mc_blas_native_srotm -

MC_TARGET_FUNC void mc_blas_native_srotm(const int n, float * x, const int incx, float * y, const int incy, const float param[5])
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	float p[5];
	mc_blas_vector_at(p, 1) = mc_blas_vector_at(param, 1);
	mc_blas_vector_at(p, 2) = mc_blas_vector_at(param, 2);
	mc_blas_vector_at(p, 3) = mc_blas_vector_at(param, 4);
	mc_blas_vector_at(p, 4) = mc_blas_vector_at(param, 3);
	mc_blas_vector_at(p, 5) = mc_blas_vector_at(param, 5);
#	else
	const float * p = param;
#	endif

#	if MC_TARGET_CPP98
	::cblas_srotm(n, x, incx, y, incy, p);
#	else
	cblas_srotm(n, x, incx, y, incy, p);
#	endif
}

#pragma mark - mc_blas_native_drotm -

MC_TARGET_FUNC void mc_blas_native_drotm(const int n, double * x, const int incx, double * y, const int incy, const double param[5])
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	double p[5];
	mc_blas_vector_at(p, 1) = mc_blas_vector_at(param, 1);
	mc_blas_vector_at(p, 2) = mc_blas_vector_at(param, 2);
	mc_blas_vector_at(p, 3) = mc_blas_vector_at(param, 4);
	mc_blas_vector_at(p, 4) = mc_blas_vector_at(param, 3);
	mc_blas_vector_at(p, 5) = mc_blas_vector_at(param, 5);
#	else
	const double * p = param;
#	endif

#	if MC_TARGET_CPP98
	::cblas_drotm(n, x, incx, y, incy, p);
#	else
	cblas_drotm(n, x, incx, y, incy, p);
#	endif
}

#endif /* !MC_BLAS_NATIVE_ROTM_H */

/* EOF */