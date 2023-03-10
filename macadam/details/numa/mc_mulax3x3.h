//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulax3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_addmulax3x3.h>

#ifndef MC_MULAX3X3_H
#define MC_MULAX3X3_H

#pragma mark - mc_mulax3x3 -

MC_TARGET_FUNC void mc_mulax3x3f(float * b, const float a[9], const float x[3])
{
//!# Requires b[3 x 1], a[3 * 3] and x[3 x 1].
//!# b=a*x
	const float x0 = x[0];
	const float x1 = x[1];
	const float x2 = x[2];

	b[0] = a[0] * x0;
	b[0] = b[0] + (a[1] * x1);
	b[0] = b[0] + (a[2] * x2);

	b[1] = a[3] * x0;
	b[1] = b[1] + (a[4] * x1);
	b[1] = b[1] + (a[5] * x2);

	b[2] = a[6] * x0;
	b[2] = b[2] + (a[7] * x1);
	b[2] = b[2] + (a[8] * x2);
}

MC_TARGET_FUNC void mc_mulax3x3ff(double * b, const float a[9], const float x[3])
{
//!# Requires b[3 x 1], a[3 * 3] and x[3 x 1].
//!# b=a*x
	b[0] = 0.0; b[1] = 0.0; b[2] = 0.0;
	mc_addmulax3x3ff(b, a, x);
}

MC_TARGET_FUNC void mc_mulax3x3(double * b, const double a[9], const double x[3])
{
//!# Requires b[3 x 1], a[3 * 3] and x[3 x 1].
//!# b=a*x
	const double x0 = x[0];
	const double x1 = x[1];
	const double x2 = x[2];

	b[0] = a[0] * x0;
	b[0] = b[0] + (a[1] * x1);
	b[0] = b[0] + (a[2] * x2);

	b[1] = a[3] * x0;
	b[1] = b[1] + (a[4] * x1);
	b[1] = b[1] + (a[5] * x2);

	b[2] = a[6] * x0;
	b[2] = b[2] + (a[7] * x1);
	b[2] = b[2] + (a[8] * x2);
}

MC_TARGET_FUNC void mc_mulax3x3l(long double * b, const long double a[9], const long double x[3])
{
//!# Requires b[3 x 1], a[3 * 3] and x[3 x 1].
//!# b=a*x
	const long double x0 = x[0];
	const long double x1 = x[1];
	const long double x2 = x[2];

	b[0] = a[0] * x0;
	b[0] = b[0] + (a[1] * x1);
	b[0] = b[0] + (a[2] * x2);

	b[1] = a[3] * x0;
	b[1] = b[1] + (a[4] * x1);
	b[1] = b[1] + (a[5] * x2);

	b[2] = a[6] * x0;
	b[2] = b[2] + (a[7] * x1);
	b[2] = b[2] + (a[8] * x2);
}

#endif /* !MC_MULAX3X3_H */

/* EOF */