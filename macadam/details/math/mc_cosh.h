//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cosh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_COSH_H
#define MC_COSH_H

#pragma mark - mc_cosh -

MC_TARGET_FUNC float mc_coshf(const float x)
{
#	if MC_TARGET_CPP98
	return ::coshf(x);
#	else
	return coshf(x);
#	endif
}

MC_TARGET_FUNC double mc_cosh(const double x)
{
#	if MC_TARGET_CPP98
	return ::cosh(x);
#	else
	return cosh(x);
#	endif
}

MC_TARGET_FUNC long double mc_coshl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::coshl(x);
#	else
	return coshl(x);
#	endif
}

#endif /* !MC_COSH_H */

/* EOF */