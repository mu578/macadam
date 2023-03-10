//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ldexp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_LDEXP_H
#define MC_LDEXP_H

#pragma mark - mc_ldexp -

MC_TARGET_FUNC float mc_ldexpf(const float x, const int n)
{
#	if MC_TARGET_CPP98
	return ::ldexpf(x, n);
#	else
	return ldexpf(x, n);
#	endif
}

MC_TARGET_FUNC double mc_ldexp(const double x, const int n)
{
#	if MC_TARGET_CPP98
	return ::ldexp(x, n);
#	else
	return ldexp(x, n);
#	endif
}

MC_TARGET_FUNC long double mc_ldexpl(const long double x, const int n)
{
#	if MC_TARGET_CPP98
	return ::ldexpl(x, n);
#	else
	return ldexpl(x, n);
#	endif
}

#endif /* !MC_LDEXP_H */

/* EOF */