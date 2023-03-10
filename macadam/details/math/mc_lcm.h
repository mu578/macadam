//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lcm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gcd.h>

#ifndef MC_LCM_H
#define MC_LCM_H

#pragma mark - mc_ilcm -

MC_TARGET_PROC int mc_ilcm(int m, int n)
{
	const int c = mc_igcd(m, n);
	return c != 0 ? ((m * n) / (c)) : 0;
}

#pragma mark - mc_llcm -

MC_TARGET_PROC long mc_llcm(long m, long n)
{
	const long c = mc_lgcd(m, n);
	return c != 0 ? ((m * n) / (c)) : 0;
}

#pragma mark - mc_lllcm -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_PROC long mc_lllcm(long long m, long long n)
{
	const long  long c = mc_llgcd(m, n);
	return c != 0 ? ((m * n) / (c)) : 0;
}
#	else
#	define mc_lllcm mc_llcm
#	endif

#endif /* !MC_LCM_H */

/* EOF */