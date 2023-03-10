//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_itrunc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ITRUNC_H
#define MC_ITRUNC_H

#pragma mark - mc_itrunc -

MC_TARGET_PROC int mc_itruncf(const float x)
{
	return mc_cast(int, x);
}

MC_TARGET_PROC int mc_itrunc(const double x)
{
	return mc_cast(int, x);
}

MC_TARGET_PROC int mc_itruncl(const long double x)
{
	return mc_cast(int, x);
}

#pragma mark - mc_ltrunc -

MC_TARGET_PROC long mc_ltruncf(const float x)
{
	return mc_cast(long, x);
}

MC_TARGET_PROC long mc_ltrunc(const double x)
{
	return mc_cast(long, x);
}

MC_TARGET_PROC long mc_ltruncl(const long double x)
{
	return mc_cast(long, x);
}

#pragma mark - mc_lltrunc -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_PROC long long mc_lltruncf(const float x)
{
	return mc_cast(long long, x);
}

MC_TARGET_PROC long long mc_lltrunc(const double x)
{
	return mc_cast(long long, x);
}

MC_TARGET_PROC long long mc_lltruncl(const long double x)
{
	return mc_cast(long long, x);
}
#	else
#	define mc_lltruncf mc_ltruncf
#	define mc_lltrunc  mc_ltrunc
#	define mc_lltruncl mc_ltruncl
#	endif

#endif /* !MC_ITRUNC_H */

/* EOF */