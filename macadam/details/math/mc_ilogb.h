//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ilogb.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ILOGB_H
#define MC_ILOGB_H

#pragma mark - mc_ilogb -

MC_TARGET_FUNC int mc_ilogbf(const float x)
{
#	if MC_TARGET_CPP98
	return ::ilogbf(x);
#	else
	return ilogbf(x);
#	endif
}

MC_TARGET_FUNC int mc_ilogb(const double x)
{
#	if MC_TARGET_CPP98
	return ::ilogb(x);
#	else
	return ilogb(x);
#	endif
}

MC_TARGET_FUNC int mc_ilogbl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::ilogbl(x);
#	else
	return ilogbl(x);
#	endif
}

#endif /* !MC_ILOGB_H */

/* EOF */