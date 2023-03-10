//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sinh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_SINH_H
#define MC_SINH_H

#pragma mark - mc_sinh -

MC_TARGET_FUNC float mc_sinhf(const float x)
{
#	if MC_TARGET_CPP98
	return ::sinhf(x);
#	else
	return sinhf(x);
#	endif
}

MC_TARGET_FUNC double mc_sinh(const double x)
{
#	if MC_TARGET_CPP98
	return ::sinh(x);
#	else
	return sinh(x);
#	endif
}

MC_TARGET_FUNC long double mc_sinhl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::sinhl(x);
#	else
	return sinhl(x);
#	endif
}

#endif /* !MC_SINH_H */

/* EOF */