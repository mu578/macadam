// mc_atanh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ATANH_H
#define MC_ATANH_H

#pragma mark - mc_atanh -

MC_TARGET_FUNC float mc_atanhf(const float x)
{
#	if MC_TARGET_CPP98
	return ::atanhf(x);
#	else
	return atanhf(x);
#	endif
}

MC_TARGET_FUNC double mc_atanh(const double x)
{
#	if MC_TARGET_CPP98
	return ::atanh(x);
#	else
	return atanh(x);
#	endif
}

MC_TARGET_FUNC long double mc_atanhl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::atanhl(x);
#	else
	return atanhl(x);
#	endif
}

#endif /* !MC_ATANH_H */

/* EOF */