//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trace2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_TRACE2X2_H
#define MC_TRACE2X2_H

#pragma mark - mc_trace2x2 -

MC_TARGET_FUNC float mc_trace2x2f(const float a[4])
{
	return a[0] + a[3];
}

MC_TARGET_FUNC double mc_trace2x2ff(const float a[4])
{
	return mc_cast(double, a[0]) + mc_cast(double, a[3]);
}

MC_TARGET_FUNC double mc_trace2x2(const double a[4])
{
	return a[0] + a[3];
}

MC_TARGET_FUNC long double mc_trace2x2l(const long double a[4])
{
	return a[0] + a[3];
}

#endif /* !MC_TRACE2X2_H */

/* EOF */