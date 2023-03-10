//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_asinh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ASINH_H
#define MC_ASINH_H

#pragma mark - mc_asinh -

MC_TARGET_FUNC float mc_asinhf(const float x)
{
#	if MC_TARGET_CPP98
	return ::asinhf(x);
#	else
	return asinhf(x);
#	endif
}

MC_TARGET_FUNC double mc_asinh(const double x)
{
#	if MC_TARGET_CPP98
	return ::asinh(x);
#	else
	return asinh(x);
#	endif
}

MC_TARGET_FUNC long double mc_asinhl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::asinhl(x);
#	else
	return asinhl(x);
#	endif
}

#endif /* !MC_ASINH_H */

/* EOF */