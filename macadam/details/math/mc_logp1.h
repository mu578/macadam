//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logp1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_LOGP1_H
#define MC_LOGP1_H

#pragma mark - mc_logp1 -

MC_TARGET_FUNC float mc_logp1f(const float x)
{
	if (x == 0.0f) {
		return x;
	}
#	if MC_TARGET_CPP98
	return ::log1pf(x);
#	else
	return log1pf(x);
#	endif
}

MC_TARGET_FUNC double mc_logp1(const double x)
{
	if (x == 0.0) {
		return x;
	}
#	if MC_TARGET_CPP98
	return ::log1p(x);
#	else
	return log1p(x);
#	endif
}

MC_TARGET_FUNC long double mc_logp1l(const long double x)
{
	if (x == 0.0L) {
		return x;
	}
#	if MC_TARGET_CPP98
	return ::log1pl(x);
#	else
	return log1pl(x);
#	endif
}

#endif /* !MC_LOGP1_H */

/* EOF */