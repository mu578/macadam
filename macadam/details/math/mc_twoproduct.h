//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_twoproduct.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fma.h>

#ifndef MC_TWOPRODUCT_H
#define MC_TWOPRODUCT_H

#pragma mark - mc_twoproduct -

MC_TARGET_FUNC void mc_twoproductf(const float a, const float b, float * x, float * y)
{
#	if MC_TARGET_HAVE_FMA
	*x = a * b;
	*y = mc_fmaf(a, b, -(*x));
#	else
//!#
//!# @note: Dekker's two-product is not a robust `fma` implementation.
//!#

//!# 2^12 + 1. ceil(FLT_MANT_DIG / 2.0) + 1.0.
	const float cs = mc_cast_expr(const float, 4096 + 1);

	float a1, a2, b1, b2, c;

	*x = a * b;
//!# split a -> a1,a2
	c  = cs * a;
	a1 = c - (c - a);
	a2 = a - a1;

//!# split b -> b1,b2
	c  = cs * b;
	b1 = c - (c - b);
	b2 = b - b1;

	*y = a2 * b2 - (*x - a1 * b1 - a2 * b1 - a1 * b2);
#	endif
}

MC_TARGET_FUNC void mc_twoproduct(const double a, const double b, double * x, double * y)
{
#	if MC_TARGET_HAVE_FMA
	*x = a * b;
	*y = mc_fma(a, b, -(*x));
#	else
//!#
//!# @note: Dekker's two-product is not a robust `fma` implementation.
//!#

//!# 2^27 + 1.
	const double cs = mc_cast_expr(const double, 134217728 + 1);

	double a1, a2, b1, b2, c;

	*x = a * b;
//!# split a -> a1,a2
	c  = cs * a;
	a1 = c - (c - a);
	a2 = a - a1;

//!# split b -> b1,b2
	c  = cs * b;
	b1 = c - (c - b);
	b2 = b - b1;

	*y = a2 * b2 - (*x - a1 * b1 - a2 * b1 - a1 * b2);
#	endif
}

MC_TARGET_FUNC void mc_twoproductl(const long double a, const long double b, long double * x, long double * y)
{
#	if MC_TARGET_HAVE_FMA
	*x = a * b;
	*y = mc_fmal(a, b, -(*x));
#	else
//!#
//!# @note: Dekker's two-product is not a robust `fma` implementation.
//!#
#	if MC_TARGET_HAVE_LONG_DOUBLE && (MC_TARGET_LONG_DOUBLE_TYPE == MC_TARGET_LONG_DOUBLE_X87)
//!# 2^32 + 1 -> ceil(LDBL_MANT_DIG / 2.0) + 1.0.
	const long double cs = mc_cast_expr(const long double, 4294967296 + 1);
#	elif MC_TARGET_HAVE_LONG_DOUBLE && (MC_TARGET_LONG_DOUBLE_TYPE == MC_TARGET_LONG_DOUBLE_ALIAS)
//!# 2^27 + 1 -> ceil(DBL_MANT_DIG / 2.0) + 1.0.
	const long double cs = mc_cast_expr(const long double, 134217728 + 1);
#	else
#	pragma message("Mantissa is too large. set @MC_TARGET_HAVE_FMA to 1.")
	const long double cs = MCK_NAN;
#	endif

	long double a1, a2, b1, b2, c;

	*x = a * b;
//!# split a -> a1,a2
	c  = cs * a;
	a1 = c - (c - a);
	a2 = a - a1;

//!# split b -> b1,b2
	c  = cs * b;
	b1 = c - (c - b);
	b2 = b - b1;

	*y = a2 * b2 - (*x - a1 * b1 - a2 * b1 - a1 * b2);
#	endif
}

#endif /* !MC_TWOPRODUCT_H */

/* EOF */