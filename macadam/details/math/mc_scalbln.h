//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_scalbln.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef SCALBLN_H
#define SCALBLN_H

#pragma mark - mc_scalbln -

MC_TARGET_FUNC float mc_scalblnf(const float x, const long y)
{
#	if MC_TARGET_CPP98
	return ::scalblnf(x, y);
#	else
	return scalblnf(x, y);
#	endif
}

MC_TARGET_FUNC double mc_scalbln(const double x, const long y)
{
#	if MC_TARGET_CPP98
	return ::scalbln(x, y);
#	else
	return scalbln(x, y);
#	endif
}

MC_TARGET_FUNC long double mc_scalblnl(const long double x, const long y)
{
#	if MC_TARGET_CPP98
	return ::scalblnl(x, y);
#	else
	return scalblnl(x, y);
#	endif
}

#endif /* !SCALBLN_H */

/* EOF */