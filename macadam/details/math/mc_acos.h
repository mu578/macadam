//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_acos.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ACOS_H
#define MC_ACOS_H

#pragma mark - mc_acos -

MC_TARGET_FUNC float mc_acosf(const float x)
{
#	if MC_TARGET_CPP98
	return ::acosf(x);
#	else
	return acosf(x);
#	endif
}

MC_TARGET_FUNC double mc_acos(const double x)
{
#	if MC_TARGET_CPP98
	return ::acos(x);
#	else
	return acos(x);
#	endif
}

MC_TARGET_FUNC long double mc_acosl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::acosl(x);
#	else
	return acosl(x);
#	endif
}

#endif /* !MC_ACOS_H */

/* EOF */