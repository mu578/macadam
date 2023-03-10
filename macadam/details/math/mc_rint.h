//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rint.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_RINT_H
#define MC_RINT_H

#pragma mark - mc_rint -

MC_TARGET_FUNC float mc_rintf(const float x)
{
#	if MC_TARGET_CPP98
	return ::rintf(x);
#	else
	return rintf(x);
#	endif
}

MC_TARGET_FUNC double mc_rint(const double x)
{
#	if MC_TARGET_CPP98
	return ::rint(x);
#	else
	return rint(x);
#	endif
}

MC_TARGET_FUNC long double mc_rintl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::rintl(x);
#	else
	return rintl(x);
#	endif
}

#endif /* !MC_RINT_H */

/* EOF */