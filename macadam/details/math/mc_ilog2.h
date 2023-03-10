//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ilog2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ILOG2_H
#define MC_ILOG2_H

#pragma mark - mc_ilog2 -

MC_TARGET_FUNC int mc_ilog2(const int x)
{
	if (x > 0 && x < MCLIMITS_IMAX) {
		const int l = mc_cast(const int, (sizeof(unsigned const int) * MCLIMITS_CBITS - 1));
		const int r = mc_cast(const int, MC_TARGET_CLZ(mc_cast(unsigned const int, x)));
		return (l - r);
	}
	return MCLIMITS_IMAX;
}

#pragma mark - mc_llog2 -

MC_TARGET_FUNC long mc_llog2(const long x)
{
	if (x > 0 && x < MCLIMITS_LMAX) {
		const long l = mc_cast(const long, (sizeof(const unsigned long) * MCLIMITS_CBITS - 1));
		const long r = mc_cast(const long, MC_TARGET_CLZL(mc_cast(const unsigned long, x)));
		return (l - r);
	}
	return MCLIMITS_LMAX;
}

#pragma mark - mc_lllog2 -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_FUNC long long mc_lllog2(const long long x)
{
	if (x > 0 && x < MCLIMITS_LLMAX) {
		const long long l = mc_cast(const long long, (sizeof(const unsigned long long) * MCLIMITS_CBITS - 1));
		const long long r = mc_cast(const long long, MC_TARGET_CLZLL(mc_cast(const unsigned long long, x)));
		return (l - r);
	}
	return MCLIMITS_LLMAX;
}
#	else
#	define mc_lllog2 mc_llog2
#	endif

#endif /* !MC_ILOG2_H */

/* EOF */