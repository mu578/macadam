//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lround.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_LROUND_H
#define MC_LROUND_H

#pragma mark - mc_lround -

MC_TARGET_FUNC long mc_lroundf(const float x)
{
#	if MC_TARGET_CPP98
	return ::lroundf(x);
#	else
	return lroundf(x);
#	endif
}

MC_TARGET_FUNC long mc_lround(const double x)
{
#	if MC_TARGET_CPP98
	return ::lround(x);
#	else
	return lround(x);
#	endif
}

MC_TARGET_FUNC long mc_lroundl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::lroundl(x);
#	else
	return lroundl(x);
#	endif
}

#endif /* !MC_LROUND_H */

/* EOF */