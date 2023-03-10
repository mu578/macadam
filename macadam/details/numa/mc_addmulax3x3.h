//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_addmulax3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ADDMULAX3X3_H
#define MC_ADDMULAX3X3_H

#pragma mark - mc_addmulax3x3 -

MC_TARGET_FUNC void mc_addmulax3x3f(float * b, const float a[9], const float x[3])
{
//!# b=b + a*x
	const float x0 = x[0];
	const float x1 = x[1];
	const float x2 = x[2];

	b[0] = b[0] + (a[0] * x0);
	b[0] = b[0] + (a[1] * x1);
	b[0] = b[0] + (a[2] * x2);

	b[1] = b[1] + (a[3] * x0);
	b[1] = b[1] + (a[4] * x1);
	b[1] = b[1] + (a[5] * x2);

	b[2] = b[2] + (a[6] * x0);
	b[2] = b[2] + (a[7] * x1);
	b[2] = b[2] + (a[8] * x2);
}

MC_TARGET_FUNC void mc_addmulax3x3ff(double * b, const float a[9], const float x[3])
{
//!# b=b + a*x
	b[0] = mc_cast(double, b[0]) + (mc_cast(double, a[0]) * mc_cast(double, x[0]));
	b[0] = mc_cast(double, b[0]) + (mc_cast(double, a[1]) * mc_cast(double, x[1]));
	b[0] = mc_cast(double, b[0]) + (mc_cast(double, a[2]) * mc_cast(double, x[2]));

	b[1] = mc_cast(double, b[1]) + (mc_cast(double, a[3]) * mc_cast(double, x[0]));
	b[1] = mc_cast(double, b[1]) + (mc_cast(double, a[4]) * mc_cast(double, x[1]));
	b[1] = mc_cast(double, b[1]) + (mc_cast(double, a[5]) * mc_cast(double, x[2]));

	b[2] = mc_cast(double, b[2]) + (mc_cast(double, a[6]) * mc_cast(double, x[0]));
	b[2] = mc_cast(double, b[2]) + (mc_cast(double, a[7]) * mc_cast(double, x[1]));
	b[2] = mc_cast(double, b[2]) + (mc_cast(double, a[8]) * mc_cast(double, x[2]));
}

MC_TARGET_FUNC void mc_addmulax3x3(double * b, const double a[9], const double x[3])
{
//!# b=b + a*x
	const double x0 = x[0];
	const double x1 = x[1];
	const double x2 = x[2];

	b[0] = b[0] + (a[0] * x0);
	b[0] = b[0] + (a[1] * x1);
	b[0] = b[0] + (a[2] * x2);

	b[1] = b[1] + (a[3] * x0);
	b[1] = b[1] + (a[4] * x1);
	b[1] = b[1] + (a[5] * x2);

	b[2] = b[2] + (a[6] * x0);
	b[2] = b[2] + (a[7] * x1);
	b[2] = b[2] + (a[8] * x2);
}

MC_TARGET_FUNC void mc_addmulax3x3l(long double * b, const long double a[9], const long double x[3])
{
//!# b=b + a*x
	const long double x0 = x[0];
	const long double x1 = x[1];
	const long double x2 = x[2];

	b[0] = b[0] + (a[0] * x0);
	b[0] = b[0] + (a[1] * x1);
	b[0] = b[0] + (a[2] * x2);

	b[1] = b[1] + (a[3] * x0);
	b[1] = b[1] + (a[4] * x1);
	b[1] = b[1] + (a[5] * x2);

	b[2] = b[2] + (a[6] * x0);
	b[2] = b[2] + (a[7] * x1);
	b[2] = b[2] + (a[8] * x2);
}

#endif /* !MC_ADDMULAX3X3_H */

/* EOF */