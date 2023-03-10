//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cbrti.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_absmag.h>

#ifndef MC_CBRTI_H
#define MC_CBRTI_H

#pragma mark - mc_cbrti -

MC_TARGET_FUNC float mc_cbrtif(const int x)
{
#	if MC_TARGET_CPP98
	return mc_iabs(x) < 0x1000001 ? (::cbrtf(mc_cast(float, x))) : mc_cast_expr(float, ::cbrt(mc_cast(const double, x)));
#	else
	return mc_iabs(x) < 0x1000001  ? (cbrtf(mc_cast(float, x)))  : mc_cast(float, cbrt(mc_cast(const double, x)));
#	endif
}

MC_TARGET_FUNC double mc_cbrti(const int x)
{
#	if MC_TARGET_CPP98
	return ::cbrt(mc_cast(double, x));
#	else
	return cbrt(mc_cast(double, x));
#	endif
}

MC_TARGET_FUNC long double mc_cbrtil(const int x)
{
#	if MC_TARGET_CPP98
	return ::cbrtl(mc_cast(long double, x));
#	else
	return cbrtl(mc_cast(long double, x));
#	endif
}

#endif /* !MC_CBRTI_H */

/* EOF */