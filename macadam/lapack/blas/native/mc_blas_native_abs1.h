//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_abs1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?abs1 computes |real(z)| + |imag(z)| of a complex number.

 *
 * \synopsis
 *    real-floating ?abs1(z)
 *    complex z
 *
 * \purpose
 *    ?abs1 computes |real(z)| + |imag(z)| of a complex number.
 *
 * \parameters
 *    [in] z - complex. The complex number argument.
 *
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#include <macadam/details/math/mc_fabs.h>

#ifndef MC_BLAS_NATIVE_ABS1_H
#define MC_BLAS_NATIVE_ABS1_H

#pragma mark - mc_blas_native_scabs1 -

MC_TARGET_FUNC float mc_blas_native_scabs1(const mc_complex_float_t z)
{
	float abs1 = 0.0f;
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		if MC_TARGET_CPP98
			abs1 = ::cblas_scabs1(&z);
#		else
			abs1 = cblas_scabs1(&z);
#		endif
#	else
	mc_unused(z);
#	endif
	return abs1;
}

#pragma mark - mc_blas_native_dcabs1 -

MC_TARGET_FUNC double mc_blas_native_dcabs1(const mc_complex_float_t z)
{
	float abs1 = 0.0f;
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		if MC_TARGET_CPP98
			abs1 = ::cblas_scabs1(&z);
#		else
			abs1 = cblas_scabs1(&z);
#		endif
#	else
	mc_unused(z);
#	endif
	return mc_cast(double, abs1);
}

#pragma mark - mc_blas_native_szabs1 -

MC_TARGET_FUNC float mc_blas_native_szabs1(const mc_complex_double_t z)
{
	double abs1 = 0.0;
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		if MC_TARGET_CPP98
			abs1 = ::cblas_dcabs1(&z);
#		else
			abs1 = cblas_dcabs1(&z);
#		endif
#	else
	mc_unused(z);
#	endif
	return mc_cast(float, abs1);
}

#pragma mark - mc_blas_native_dzabs1 -

MC_TARGET_FUNC double mc_blas_native_dzabs1(const mc_complex_double_t z)
{
	double abs1 = 0.0;
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		if MC_TARGET_CPP98
			abs1 = ::cblas_dcabs1(&z);
#		else
			abs1 = cblas_dcabs1(&z);
#		endif
#	else
	mc_unused(z);
#	endif
	return abs1;
}

#endif /* !MC_BLAS_NATIVE_ABS1_H */

/* EOF */