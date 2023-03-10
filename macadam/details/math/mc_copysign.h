//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_copysign.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_COPYSIGN_H
#define MC_COPYSIGN_H

#pragma mark - mc_copysign -

MC_TARGET_FUNC float mc_copysignf(const float x, const float y)
{
#	if MC_TARGET_CPP98
	return ::copysignf(x, y);
#	else
	return copysignf(x, y);
#	endif
}

MC_TARGET_FUNC double mc_copysign(const double x, const double y)
{
#	if MC_TARGET_CPP98
	return ::copysign(x, y);
#	else
	return copysign(x, y);
#	endif
}

MC_TARGET_FUNC long double mc_copysignl(const long double x, const long double y)
{
#	if MC_TARGET_CPP98
	return ::copysignl(x, y);
#	else
	return copysignl(x, y);
#	endif
}

#endif /* !MC_COPYSIGN_H */

/* EOF */