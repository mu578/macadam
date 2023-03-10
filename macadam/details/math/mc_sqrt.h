//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sqrt.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_SQRT_H
#define MC_SQRT_H

#pragma mark - mc_sqrt -

MC_TARGET_FUNC float mc_sqrtf(const float x)
{
#	if MC_TARGET_CPP98
	return ::sqrtf(x);
#	else
	return sqrtf(x);
#	endif
}

MC_TARGET_FUNC double mc_sqrt(const double x)
{
#	if MC_TARGET_CPP98
	return ::sqrt(x);
#	else
	return sqrt(x);
#	endif
}

MC_TARGET_FUNC long double mc_sqrtl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::sqrtl(x);
#	else
	return sqrtl(x);
#	endif
}

#endif /* !MC_SQRT_H */

/* EOF */