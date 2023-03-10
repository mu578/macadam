//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cbrtu.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_CBRTU_H
#define MC_CBRTU_H

#pragma mark - mc_cbrtu -

MC_TARGET_FUNC float mc_cbrtuf(const unsigned int x)
{
#	if MC_TARGET_CPP98
	return x < 0x1000001U ? (::cbrtf(mc_cast(float, x))) : mc_cast_expr(float, ::cbrt(mc_cast(double, x)));
#	else
	return x < 0x1000001U  ? (cbrtf(mc_cast(float, x)))  : mc_cast(float, cbrt(mc_cast(double, x)));
#	endif
}

MC_TARGET_FUNC double mc_cbrtu(const unsigned int x)
{
#	if MC_TARGET_CPP98
	return ::cbrt(mc_cast(double, x));
#	else
	return cbrt(mc_cast(double, x));
#	endif
}

MC_TARGET_FUNC long double mc_cbrtul(const unsigned int x)
{
#	if MC_TARGET_CPP98
	return ::cbrtl(mc_cast(long double, x));
#	else
	return cbrtl(mc_cast(long double, x));
#	endif
}

#endif /* !MC_CBRTU_H */

/* EOF */