//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_det2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_DET2X2_H
#define MC_DET2X2_H

#pragma mark - mc_det2x2 -

MC_TARGET_FUNC float mc_det2x2f(const float a[4])
{
	return a[0] * a[3] - a[2] * a[1];
}

MC_TARGET_FUNC double mc_det2x2ff(const float a[4])
{
	return mc_cast(double, a[0]) * mc_cast(double, a[3]) - mc_cast(double, a[2]) * mc_cast(double, a[1]);
}

MC_TARGET_FUNC double mc_det2x2(const double a[4])
{
	return a[0] * a[3] - a[2] * a[1];
}

MC_TARGET_FUNC long double mc_det2x2l(const long double a[4])
{
	return a[0] * a[3] - a[2] * a[1];
}

#endif /* !MC_DET2X2_H */

/* EOF */