//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lrint.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_LRINT_H
#define MC_LRINT_H

#pragma mark - mc_lrint -

MC_TARGET_FUNC long mc_lrintf(const float x)
{
#	if MC_TARGET_CPP98
	return ::lrintf(x);
#	else
	return lrintf(x);
#	endif
}

MC_TARGET_FUNC long mc_lrint(const double x)
{
#	if MC_TARGET_CPP98
	return ::lrint(x);
#	else
	return lrint(x);
#	endif
}

MC_TARGET_FUNC long mc_lrintl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::lrintl(x);
#	else
	return lrintl(x);
#	endif
}

#endif /* !MC_LRINT_H */

/* EOF */