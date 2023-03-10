//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ceil.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_CEIL_H
#define MC_CEIL_H

#pragma mark - mc_ceil -

MC_TARGET_FUNC float mc_ceilf(const float x)
{
#	if MC_TARGET_CPP98
	return ::ceilf(x);
#	else
	return ceilf(x);
#	endif
}

MC_TARGET_FUNC double mc_ceil(const double x)
{
#	if MC_TARGET_CPP98
	return ::ceil(x);
#	else
	return ceil(x);
#	endif
}

MC_TARGET_FUNC long double mc_ceill(const long double x)
{
#	if MC_TARGET_CPP98
	return ::ceill(x);
#	else
	return ceill(x);
#	endif
}

#endif /* !MC_CEIL_H */

/* EOF */