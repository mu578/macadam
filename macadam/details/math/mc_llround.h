//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_llround.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_lround.h>

#ifndef MC_LLROUND_H
#define MC_LLROUND_H

#pragma mark - mc_llround -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_FUNC long long mc_llroundf(const float x)
{
#	if MC_TARGET_CPP11
	return ::llroundf(x);
#	else
	return llroundf(x);
#	endif
}

MC_TARGET_FUNC long long mc_llround(const double x)
{
#	if MC_TARGET_CPP11
	return ::llround(x);
#	else
	return llround(x);
#	endif
}

MC_TARGET_FUNC long long mc_llroundl(const long double x)
{
#	if MC_TARGET_CPP11
	return ::llroundl(x);
#	else
	return llroundl(x);
#	endif
}
#	else
#	define mc_llroundf mc_lroundf
#	define mc_llround  mc_lround
#	define mc_llroundl mc_lroundl
#	endif

#endif /* !MC_LLROUND_H */

/* EOF */