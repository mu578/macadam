// mc_tan.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_TAN_H
#define MC_TAN_H

#pragma mark - mc_tan -

MC_TARGET_FUNC float mc_tanf(const float x)
{
#	if MC_TARGET_CPP98
	return ::tanf(x);
#	else
	return tanf(x);
#	endif
}

MC_TARGET_FUNC double mc_tan(const double x)
{
#	if MC_TARGET_CPP98
	return ::tan(x);
#	else
	return tan(x);
#	endif
}

MC_TARGET_FUNC long double mc_tanl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::tanl(x);
#	else
	return tanl(x);
#	endif
}

#endif /* !MC_TAN_H */

/* EOF */