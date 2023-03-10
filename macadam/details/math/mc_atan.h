//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_atan.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ATAN_H
#define MC_ATAN_H

#pragma mark - mc_atan -

MC_TARGET_FUNC float mc_atanf(const float x)
{
#	if MC_TARGET_CPP98
	return ::atanf(x);
#	else
	return atanf(x);
#	endif
}

MC_TARGET_FUNC double mc_atan(const double x)
{
#	if MC_TARGET_CPP98
	return ::atan(x);
#	else
	return atan(x);
#	endif
}

MC_TARGET_FUNC long double mc_atanl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::atanl(x);
#	else
	return atanl(x);
#	endif
}

#endif /* !MC_ATAN_H */

/* EOF */