//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_nchoosek.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_NCHOOSEK_H
#define MC_NCHOOSEK_H

#pragma mark - mc_nchoosek -

MC_TARGET_FUNC unsigned int mc_nchoosek(const unsigned int n, const unsigned int k)
{
	if (n < MCLIMITS_UIMAX && k < MCLIMITS_UIMAX) {
		if (k > n) {
			return MCLIMITS_UIMAX;
		} else if ((k == 0) || (k == n)) {
			return 1;
		} else if ((k == 1) || (k == n - 1)) {
			return n;
		}
		return (n * mc_nchoosek(n - 1, k - 1)) / k;
	}
	return MCLIMITS_UIMAX;
}

#pragma mark - mc_nchoosekul -

MC_TARGET_PROC unsigned long mc_nchoosekul(const unsigned long n, const unsigned long k)
{
	if (n < MCLIMITS_ULMAX && k < MCLIMITS_ULMAX) {
		if (k > n) {
			return MCLIMITS_ULMAX;
		} else if ((k == 0) || (k == n)) {
			return 1;
		} else if ((k == 1) || (k == n - 1)) {
			return n;
		}
		return (n * mc_nchoosekul(n - 1, k - 1)) / k;
	}
	return MCLIMITS_ULMAX;
}

#pragma mark - mc_nchoosekull -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_PROC unsigned long long mc_nchoosekull(const unsigned long long n, const unsigned long long k)
{
	if (n < MCLIMITS_ULLMAX && k < MCLIMITS_ULLMAX) {
		if (k > n) {
			return MCLIMITS_ULLMAX;
		} else if ((k == 0) || (k == n)) {
			return 1;
		} else if ((k == 1) || (k == n - 1)) {
			return n;
		}
		return (n * mc_nchoosekull(n - 1, k - 1)) / k;
	}
	return MCLIMITS_ULLMAX;
}
#	else
#	define mc_nchoosekull mc_nchoosekul
#	endif

#endif /* !MC_NCHOOSEK_H */

/* EOF */