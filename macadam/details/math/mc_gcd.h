//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_gcd.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_GCD_H
#define MC_GCD_H

#pragma mark - mc_igcd -

MC_TARGET_PROC int mc_igcd(const int m, const int n)
{
	if (m == 0 && n == 0) {
		return 0;
	}
	if (n == 0) {
		return m;
	}
	int a = m < 0 ? -m : m;
	int b = n < 0 ? -n : n;
	int c = a % b;
	while (c > 0) { a = b; b = c; c = a % b; }
	return b;
}

#pragma mark - mc_lgcd -

MC_TARGET_PROC long mc_lgcd(const long m, const long n)
{
	if (m == 0 && n == 0) {
		return 0;
	}
	if (n == 0) {
		return m;
	}
	long a = m < 0 ? -m : m;
	long b = n < 0 ? -n : n;
	long c = a % b;
	while (c > 0) { a = b; b = c; c = a % b; }
	return b;
}

#pragma mark - mc_llgcd -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_PROC long long mc_llgcd(const long long m, const long long n)
{
	if (m == 0 && n == 0) {
		return 0;
	}
	if (n == 0) {
		return m;
	}
	long long a = m < 0 ? -m : m;
	long long b = n < 0 ? -n : n;
	long long c = a % b;
	while (c > 0) { a = b; b = c; c = a % b; }
	return b;
}
#	else
#	define mc_llgcd mc_lgcd
#	endif

#endif /* !MC_GCD_H */

/* EOF */