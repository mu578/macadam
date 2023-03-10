//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_det3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_DET3X3_H
#define MC_DET3X3_H

#pragma mark - mc_det3x3 -

MC_TARGET_FUNC float mc_det3x3f(const float a[9])
{
	return a[0] * ((a[4] * a[8]) - (a[7] * a[5]))
		  - a[1] * ( a[3] * a[8]  -  a[6] * a[5] )
		  + a[2] * ( a[3] * a[7]  -  a[6] * a[4] );
}

MC_TARGET_FUNC double mc_det3x3ff(const float a[9])
{
	return mc_cast(double, a[0]) * ((mc_cast(double, a[4]) * mc_cast(double, a[8])) - (mc_cast(double, a[7]) * mc_cast(double, a[5])))
		  - mc_cast(double, a[1]) * ( mc_cast(double, a[3]) * mc_cast(double, a[8])  -  mc_cast(double, a[6]) * mc_cast(double, a[5]) )
		  + mc_cast(double, a[2]) * ( mc_cast(double, a[3]) * mc_cast(double, a[7])  -  mc_cast(double, a[6]) * mc_cast(double, a[4]) );
}

MC_TARGET_FUNC double mc_det3x3(const double a[9])
{
	return a[0] * ((a[4] * a[8]) - (a[7] * a[5]))
		  - a[1] * ( a[3] * a[8]  -  a[6] * a[5] )
		  + a[2] * ( a[3] * a[7]  -  a[6] * a[4] );
}

MC_TARGET_FUNC long double mc_det3x3l(const long double a[9])
{
	return a[0] * ((a[4] * a[8]) - (a[7] * a[5]))
		  - a[1] * ( a[3] * a[8]  -  a[6] * a[5] )
		  + a[2] * ( a[3] * a[7]  -  a[6] * a[4] );
}

#endif /* !MC_DET3X3_H */

/* EOF */