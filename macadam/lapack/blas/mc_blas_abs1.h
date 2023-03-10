//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_abs1.h
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
 *    [in] z - complex. The complex argument.
 *
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#include <macadam/details/math/mc_fabs.h>

#ifndef MC_BLAS_ABS1_H
#define MC_BLAS_ABS1_H

#pragma mark - mc_blas_scabs1 -

MC_TARGET_FUNC float mc_blas_scabs1(const mc_complex_float_t z)
{
	return mc_fabsf(mc_cmplxrf(z)) + mc_fabsf(mc_cmplxif(z));
}

#pragma mark - mc_blas_dcabs1 -

MC_TARGET_FUNC double mc_blas_dcabs1(const mc_complex_float_t z)
{
	return mc_cast_expr(double, mc_fabsf(mc_cmplxrf(z)) + mc_fabsf(mc_cmplxif(z)));
}

#pragma mark - mc_blas_lcabs1 -

MC_TARGET_FUNC long double mc_blas_lcabs1(const mc_complex_float_t z)
{
	return mc_cast_expr(long double, mc_fabsf(mc_cmplxrf(z)) + mc_fabsf(mc_cmplxif(z)));
}

#pragma mark - mc_blas_szabs1 -

MC_TARGET_FUNC float mc_blas_szabs1(const mc_complex_double_t z)
{
	return mc_cast_expr(float, mc_fabs(mc_cmplxr(z)) + mc_fabs(mc_cmplxi(z)));
}

#pragma mark - mc_blas_dzabs1 -

MC_TARGET_FUNC double mc_blas_dzabs1(const mc_complex_double_t z)
{
	return mc_fabs(mc_cmplxr(z)) + mc_fabs(mc_cmplxi(z));
}

#pragma mark - mc_blas_lzabs1 -

MC_TARGET_FUNC long double mc_blas_lzabs1(const mc_complex_double_t z)
{
	return mc_cast_expr(long double, mc_fabs(mc_cmplxr(z)) + mc_fabs(mc_cmplxi(z)));
}

#pragma mark - mc_blas_sqabs1 -

MC_TARGET_FUNC float mc_blas_sqabs1(const mc_complex_long_double_t z)
{
	return mc_cast_expr(float, mc_fabsl(mc_cmplxrl(z)) + mc_fabsl(mc_cmplxil(z)));
}

#pragma mark - mc_blas_dqabs1 -

MC_TARGET_FUNC double mc_blas_dqabs1(const mc_complex_long_double_t z)
{
	return mc_cast_expr(double, mc_fabsl(mc_cmplxrl(z)) + mc_fabsl(mc_cmplxil(z)));
}

#pragma mark - mc_blas_lqabs1 -

MC_TARGET_FUNC long double mc_blas_lqabs1(const mc_complex_long_double_t z)
{
	return mc_fabsl(mc_cmplxrl(z)) + mc_fabsl(mc_cmplxil(z));
}

#endif /* !MC_BLAS_ABS1_H */

/* EOF */