//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fdim.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_FDIM_H
#define MC_FDIM_H

#pragma mark - mc_fdim -

MC_TARGET_FUNC float mc_fdimf(const float x, const float y)
{
#	if MC_TARGET_CPP98
	return ::fdimf(x, y);
#	else
	return fdimf(x, y);
#	endif
}

MC_TARGET_FUNC double mc_fdim(const double x, const double y)
{
#	if MC_TARGET_CPP98
	return ::fdim(x, y);
#	else
	return fdim(x, y);
#	endif
}

MC_TARGET_FUNC long double mc_fdiml(const long double x, const long double y)
{
#	if MC_TARGET_CPP98
	return ::fdiml(x, y);
#	else
	return fdiml(x, y);
#	endif
}

#endif /* !MC_FDIM_H */

/* EOF */