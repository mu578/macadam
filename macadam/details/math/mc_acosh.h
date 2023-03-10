//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_acosh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ACOSH_H
#define MC_ACOSH_H

#pragma mark - mc_acosh -

MC_TARGET_FUNC float mc_acoshf(const float x)
{
#	if MC_TARGET_CPP98
	return ::acoshf(x);
#	else
	return acoshf(x);
#	endif
}

MC_TARGET_FUNC double mc_acosh(const double x)
{
#	if MC_TARGET_CPP98
	return ::acosh(x);
#	else
	return acosh(x);
#	endif
}

MC_TARGET_FUNC long double mc_acoshl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::acoshl(x);
#	else
	return acoshl(x);
#	endif
}

#endif /* !MC_ACOSH_H */

/* EOF */