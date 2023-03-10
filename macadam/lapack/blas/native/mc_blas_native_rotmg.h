//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_rotmg.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?rotmg constructs a Gentleman's modified Given's plane rotation.
 *
 * \synopsis
 *    real-floating ?rotmg(d1, d2, x1, y1, param)
 *    real-floating d1, d2, x1, y1
 *    real-floating param(5)
 *
 * \purpose
 *    ?rotmg computes the parameters for a modified Givens rotation. Construct Gentleman's
 *    modified a Given's plane rotation that will annihilate an element of a vector:
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
 * \parameters
 *    [out] d1    - real-floating. The first diagonal entry in the H matrix. Overwritten to
 *    reflect the effect of the transformation.
 *
 *    [out] d2    - real-floating. The second diagonal entry in the H matrix. Overwritten to
 *    reflect the effect of the transformation.
 *
 *    [out] x1    - real-floating. The first element of the vector to which the H matrix is
 *    applied.  to reflect the effect of the transformation.
 *
 *    [in] y1     - real-floating. The second element of the vector to which the H matrix is
 *    applied.
 *
 *    [out] param - real-floating array of size 5 such as param[flag, h11, h21, h12, h22].
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

#ifndef MC_BLAS_NATIVE_ROTMG_H
#define MC_BLAS_NATIVE_ROTMG_H

#pragma mark - mc_blas_native_srotmg -

MC_TARGET_FUNC void mc_blas_native_srotmg(float * d1, float * d2, float * x1, const float y1, float param[5])
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	float w;
#	endif

#	if MC_TARGET_CPP98
	::cblas_srotmg(d1, d2, x1, y1, param);
#	else
	cblas_srotmg(d1, d2, x1, y1, param);
#	endif

#	if MC_TARGET_BLAS_USE_CLAYOUT
	mcswap_var(w,
		  mc_blas_vector_at(param, 3)
		, mc_blas_vector_at(param, 4)
	);
#	endif
}

#pragma mark - mc_blas_native_drotmg -

MC_TARGET_FUNC void mc_blas_native_drotmg(double * d1, double * d2, double * x1, const double y1, double param[5])
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	double w;
#	endif

#	if MC_TARGET_CPP98
	::cblas_drotmg(d1, d2, x1, y1, param);
#	else
	cblas_drotmg(d1, d2, x1, y1, param);
#	endif

#	if MC_TARGET_BLAS_USE_CLAYOUT
	mcswap_var(w
		, mc_blas_vector_at(param, 3)
		, mc_blas_vector_at(param, 4)
	);
#	endif
}

#endif /* !MC_BLAS_NATIVE_ROTMG_H */

/* EOF */