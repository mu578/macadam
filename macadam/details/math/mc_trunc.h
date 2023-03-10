//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trunc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_TRUNC_H
#define MC_TRUNC_H

#pragma mark - mc_trunc -

MC_TARGET_FUNC float mc_truncf(const float x)
{
#	if MC_TARGET_CPP98
	return ::truncf(x);
#	else
	return truncf(x);
#	endif
}

MC_TARGET_FUNC double mc_trunc(const double x)
{
#	if MC_TARGET_CPP98
	return ::trunc(x);
#	else
	return trunc(x);
#	endif
}

MC_TARGET_FUNC long double mc_truncl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::truncl(x);
#	else
	return truncl(x);
#	endif
}

#endif /* !MC_TRUNC_H */

/* EOF */