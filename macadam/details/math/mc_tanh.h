// mc_tanh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_TANH_H
#define MC_TANH_H

#pragma mark - mc_tanh -

MC_TARGET_FUNC float mc_tanhf(const float x)
{
#	if MC_TARGET_CPP98
	return ::tanhf(x);
#	else
	return tanhf(x);
#	endif
}

MC_TARGET_FUNC double mc_tanh(const double x)
{
#	if MC_TARGET_CPP98
	return ::tanh(x);
#	else
	return tanh(x);
#	endif
}

MC_TARGET_FUNC long double mc_tanhl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::tanhl(x);
#	else
	return tanhl(x);
#	endif
}

#endif /* !MC_TANH_H */

/* EOF */