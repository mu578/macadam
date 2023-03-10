//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trace3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_TRACE3X3_H
#define MC_TRACE3X3_H

#pragma mark - mc_trace3x3 -

MC_TARGET_FUNC float mc_trace3x3f(const float a[9])
{
	return a[0] + a[4] + a[8];
}

MC_TARGET_FUNC double mc_trace3x3ff(const float a[9])
{
	return mc_cast(double, a[0]) + mc_cast(double, a[4]) + mc_cast(double, a[8]);
}

MC_TARGET_FUNC double mc_trace3x3(const double a[9])
{
	return a[0] + a[4] + a[8];
}

MC_TARGET_FUNC long double mc_trace3x3l(const long double a[9])
{
	return a[0] + a[4] + a[8];
}

#endif /* !MC_TRACE3X3_H */

/* EOF */