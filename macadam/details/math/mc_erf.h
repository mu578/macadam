//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_erf.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ERF_H
#define MC_ERF_H

#pragma mark - mc_erf -

MC_TARGET_FUNC float mc_erff(const float x)
{
#	if MC_TARGET_CPP98
	return ::erff(x);
#	else
	return erff(x);
#	endif
}

MC_TARGET_FUNC double mc_erf(const double x)
{
#	if MC_TARGET_CPP98
	return ::erf(x);
#	else
	return erf(x);
#	endif
}

MC_TARGET_FUNC long double mc_erfl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::erfl(x);
#	else
	return erfl(x);
#	endif
}

#endif /* !MC_ERF_H */

/* EOF */