//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sin.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_SIN_H
#define MC_SIN_H

#pragma mark - mc_sin -

MC_TARGET_FUNC float mc_sinf(const float x)
{
#	if MC_TARGET_CPP98
	return ::sinf(x);
#	else
	return sinf(x);
#	endif
}

MC_TARGET_FUNC double mc_sin(const double x)
{
#	if MC_TARGET_CPP98
	return ::sin(x);
#	else
	return sin(x);
#	endif
}

MC_TARGET_FUNC long double mc_sinl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::sinl(x);
#	else
	return sinl(x);
#	endif
}

#endif /* !MC_SIN_H */

/* EOF */