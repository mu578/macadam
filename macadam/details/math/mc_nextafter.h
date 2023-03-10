//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_nextafter.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_NEXTAFTER_H
#define MC_NEXTAFTER_H

#pragma mark - mc_nextafter -

MC_TARGET_FUNC float mc_nextafterf(const float x, const float y)
{
#	if MC_TARGET_CPP98
	return ::nextafterf(x, y);
#	else
	return nextafterf(x, y);
#	endif
}

MC_TARGET_FUNC double mc_nextafter(const double x, const double y)
{
#	if MC_TARGET_CPP98
	return ::nextafter(x, y);
#	else
	return nextafter(x, y);
#	endif
}

MC_TARGET_FUNC long double mc_nextafterl(const long double x, const long double y)
{
#	if MC_TARGET_CPP98
	return ::nextafterl(x, y);
#	else
	return nextafterl(x, y);
#	endif
}

#endif /* !MC_NEXTAFTER_H */

/* EOF */