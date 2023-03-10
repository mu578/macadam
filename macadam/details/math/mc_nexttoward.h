//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_nexttoward.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_NEXTTOWARD_H
#define MC_NEXTTOWARD_H

#pragma mark - mc_nexttoward -

MC_TARGET_FUNC float mc_nexttowardf(const float x, const float y)
{
#	if MC_TARGET_CPP98
	return ::nexttowardf(x, y);
#	else
	return nexttowardf(x, y);
#	endif
}

MC_TARGET_FUNC double mc_nexttoward(const double x, const double y)
{
#	if MC_TARGET_CPP98
	return ::nexttoward(x, y);
#	else
	return nexttoward(x, y);
#	endif
}

MC_TARGET_FUNC long double mc_nexttowardl(const long double x, const long double y)
{
#	if MC_TARGET_CPP98
	return ::nexttowardl(x, y);
#	else
	return nexttowardl(x, y);
#	endif
}

#endif /* !MC_NEXTTOWARD_H */

/* EOF */