//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulab2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MULAB2X2_H
#define MC_MULAB2X2_H

#pragma mark - mc_mulab2x2 -

MC_TARGET_FUNC void mc_mulab2x2f(float * c, const float a[4], const float b[4])
{
//!# c=a*b
	c[0] = (a[0] * b[0]) + (a[1] * b[2]);
	c[1] = (a[0] * b[1]) + (a[1] * b[3]);

	c[2] = (a[2] * b[0]) + (a[3] * b[2]);
	c[3] = (a[2] * b[1]) + (a[3] * b[3]);
}

MC_TARGET_FUNC void mc_mulab2x2ff(double * c, const float a[4], const float b[4])
{
//!# c=a*b
	c[0] = (mc_cast(double, a[0]) * mc_cast(double, b[0])) + (mc_cast(double, a[1]) * mc_cast(double, b[2]));
	c[1] = (mc_cast(double, a[0]) * mc_cast(double, b[1])) + (mc_cast(double, a[1]) * mc_cast(double, b[3]));

	c[2] = (mc_cast(double, a[2]) * mc_cast(double, b[0])) + (mc_cast(double, a[3]) * mc_cast(double, b[2]));
	c[3] = (mc_cast(double, a[2]) * mc_cast(double, b[1])) + (mc_cast(double, a[3]) * mc_cast(double, b[3]));
}

MC_TARGET_FUNC void mc_mulab2x2fd(double * MC_TARGET_RESTRICT c, const float a[4], const double b[4])
{
//!# c=a*b
	c[0] = (mc_cast(double, a[0]) * b[0]) + (mc_cast(double, a[1]) * b[2]);
	c[1] = (mc_cast(double, a[0]) * b[1]) + (mc_cast(double, a[1]) * b[3]);

	c[2] = (mc_cast(double, a[2]) * b[0]) + (mc_cast(double, a[3]) * b[2]);
	c[3] = (mc_cast(double, a[2]) * b[1]) + (mc_cast(double, a[3]) * b[3]);
}

MC_TARGET_FUNC void mc_mulab2x2df(double * MC_TARGET_RESTRICT c, const double a[4], const float b[4])
{
//!# c=a*b
	c[0] = (a[0] * mc_cast(double, b[0])) + (a[1] * mc_cast(double, b[2]));
	c[1] = (a[0] * mc_cast(double, b[1])) + (a[1] * mc_cast(double, b[3]));

	c[2] = (a[2] * mc_cast(double, b[0])) + (a[3] * mc_cast(double, b[2]));
	c[3] = (a[2] * mc_cast(double, b[1])) + (a[3] * mc_cast(double, b[3]));
}

MC_TARGET_FUNC void mc_mulab2x2(double * c, const double a[4], const double b[4])
{
//!# c=a*b
	c[0] = (a[0] * b[0]) + (a[1] * b[2]);
	c[1] = (a[0] * b[1]) + (a[1] * b[3]);

	c[2] = (a[2] * b[0]) + (a[3] * b[2]);
	c[3] = (a[2] * b[1]) + (a[3] * b[3]);
}

MC_TARGET_FUNC void mc_mulab2x2l(long double * c, const long double a[4], const long double b[4])
{
//!# c=a*b
	c[0] = (a[0] * b[0]) + (a[1] * b[2]);
	c[1] = (a[0] * b[1]) + (a[1] * b[3]);

	c[2] = (a[2] * b[0]) + (a[3] * b[2]);
	c[3] = (a[2] * b[1]) + (a[3] * b[3]);
}

#endif /* !MC_MULAB2X2_H */

/* EOF */