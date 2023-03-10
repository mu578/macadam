//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_atan2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ATAN2_H
#define MC_ATAN2_H

#pragma mark - mc_atan2 -

MC_TARGET_FUNC float mc_atan2f(const float y, const float x)
{
#	if MC_TARGET_CPP98
	return ::atan2f(y, x);
#	else
	return atan2f(y, x);
#	endif
}

MC_TARGET_FUNC double mc_atan2(const double y, const double x)
{
#	if MC_TARGET_CPP98
	return ::atan2(y, x);
#	else
	return atan2(y, x);
#	endif
}

MC_TARGET_FUNC long double mc_atan2l(const long double y, const long double x)
{
#	if MC_TARGET_CPP98
	return ::atan2l(y, x);
#	else
	return atan2l(y, x);
#	endif
}

#endif /* !MC_ATAN2_H */

/* EOF */