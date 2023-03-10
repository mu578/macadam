//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_iexp2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_IEXP2_H
#define MC_IEXP2_H

#pragma mark - mc_iexp2 -

MC_TARGET_FUNC int mc_iexp2(const int e)
{
	const int argmax = sizeof(int);
	const int retmax = MCLIMITS_IMAX;
	return !(e < 0 || e >= argmax) ? 1 << e : retmax;
}

#pragma mark - mc_lexp2 -

MC_TARGET_FUNC long mc_lexp2(const long e)
{
	const long argmax = sizeof(long);
	const long retmax = MCLIMITS_LMAX;
	return !(e < 0 || e >= argmax) ? 1 << e : retmax;
}

#pragma mark - mc_llexp2 -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_FUNC long long mc_llexp2(const long long e)
{
	const long long argmax = sizeof(long long);
	const long long retmax = MCLIMITS_LLMAX;
	return !(e < 0 || e >= argmax) ? 1 << e : retmax;
}
#	else
#	define mc_llexp2 mc_lexp2
#	endif

#endif /* !MC_IEXP2_H */

/* EOF */