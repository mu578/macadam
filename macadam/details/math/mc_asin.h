//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_asin.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ASIN_H
#define MC_ASIN_H

#pragma mark - mc_asin -

MC_TARGET_FUNC float mc_asinf(const float x)
{
#	if MC_TARGET_CPP98
	return ::asinf(x);
#	else
	return asinf(x);
#	endif
}

MC_TARGET_FUNC double mc_asin(const double x)
{
#	if MC_TARGET_CPP98
	return ::asin(x);
#	else
	return asin(x);
#	endif
}

MC_TARGET_FUNC long double mc_asinl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::asinl(x);
#	else
	return asinl(x);
#	endif
}

#endif /* !MC_ASIN_H */

/* EOF */