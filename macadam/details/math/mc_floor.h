//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_floor.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_FLOOR_H
#define MC_FLOOR_H

#pragma mark - mc_floor -

MC_TARGET_FUNC float mc_floorf(const float x)
{
#	if MC_TARGET_CPP98
	return ::floorf(x);
#	else
	return floorf(x);
#	endif
}

MC_TARGET_FUNC double mc_floor(const double x)
{
#	if MC_TARGET_CPP98
	return ::floor(x);
#	else
	return floor(x);
#	endif
}

MC_TARGET_FUNC long double mc_floorl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::floorl(x);
#	else
	return floorl(x);
#	endif
}

#endif /* !MC_FLOOR_H */

/* EOF */