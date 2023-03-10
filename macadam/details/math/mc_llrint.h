//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_llrint.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_lrint.h>

#ifndef MC_LLRINT_H
#define MC_LLRINT_H

#pragma mark - mc_llrint -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_FUNC long long mc_llrintf(const float x)
{
#	if MC_TARGET_CPP11
	return ::llrintf(x);
#	else
	return llrintf(x);
#	endif
}

MC_TARGET_FUNC long long mc_llrint(const double x)
{
#	if MC_TARGET_CPP11
	return ::llrint(x);
#	else
	return llrint(x);
#	endif
}

MC_TARGET_FUNC long long mc_llrintl(const long double x)
{
#	if MC_TARGET_CPP11
	return ::llrintl(x);
#	else
	return llrintl(x);
#	endif
}
#	else
#	define mc_llrintf mc_lrintf
#	define mc_llrint  mc_lrint
#	define mc_llrintl mc_lrintl
#	endif

#endif /* !MC_LLRINT_H */

/* EOF */