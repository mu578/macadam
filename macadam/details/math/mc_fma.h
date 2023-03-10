//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_FMA_H
#define MC_FMA_H

#pragma mark - mc_fma -

MC_TARGET_FUNC float mc_fmaf(const float x, const float y, const float z)
{
#	if MC_TARGET_CPP98
	return ::fmaf(x, y, z);
#	else
	return fmaf(x, y, z);
#	endif
}

MC_TARGET_FUNC double mc_fmaff(const float x, const float y, const float z)
{
#	if MC_TARGET_CPP98
	return ::fma(mc_cast(const double, x), mc_cast(const double, y), mc_cast(const double, z));
#	else
	return fma(mc_cast(const double, x), mc_cast(const double, y), mc_cast(const double, z));
#	endif
}

MC_TARGET_FUNC double mc_fmafd(const float x, const float y, const double z)
{
#	if MC_TARGET_CPP98
	return ::fma(mc_cast(const double, x), mc_cast(const double, y), z);
#	else
	return fma(mc_cast(const double, x), mc_cast(const double, y), z);
#	endif
}

MC_TARGET_FUNC double mc_fmadf(const double x, const double y, const float z)
{
#	if MC_TARGET_CPP98
	return ::fma(x, y, mc_cast(const double, z));
#	else
	return fma(x, y, mc_cast(const double, z));
#	endif
}

MC_TARGET_FUNC double mc_fma(const double x, const double y, const double z)
{
#	if MC_TARGET_CPP98
	return ::fma(x, y, z);
#	else
	return fma(x, y, z);
#	endif
}

MC_TARGET_FUNC long double mc_fmal(const long double x, const long double y, const long double z)
{
#	if MC_TARGET_CPP98
	return ::fmal(x, y, z);
#	else
	return fmal(x, y, z);
#	endif
}

#endif /* !MC_FMA_H */

/* EOF */