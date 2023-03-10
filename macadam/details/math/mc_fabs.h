//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fabs.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_FABS_H
#define MC_FABS_H

#pragma mark - mc_fabs -

MC_TARGET_FUNC float mc_fabsf(const float x)
{
#	if MC_TARGET_CPP98
	return ::fabsf(x);
#	else
	return fabsf(x);
#	endif
}

MC_TARGET_FUNC double mc_fabs(const double x)
{
#	if MC_TARGET_CPP98
	return ::fabs(x);
#	else
	return fabs(x);
#	endif
}

MC_TARGET_FUNC long double mc_fabsl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::fabsl(x);
#	else
	return fabsl(x);
#	endif
}

#endif /* !MC_FABS_H */

/* EOF */