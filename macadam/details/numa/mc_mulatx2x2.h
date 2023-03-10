//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulatx2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_addmulatx2x2.h>

#ifndef MC_MULATX2X2_H
#define MC_MULATX2X2_H

#pragma mark - mc_mulatx2x2 -

MC_TARGET_FUNC void mc_mulatx2x2f(float * b, const float a[4], const float x[2])
{
//!# b=a'*x
	const float x0 = x[0];
	const float x1 = x[1];

	b[0] = a[0] * x0;
	b[0] = b[0] + (a[2] * x1);

	b[1] = a[1] * x0;
	b[1] = b[1] + (a[3] * x1);
}

MC_TARGET_FUNC void mc_mulatx2x2ff(double * b, const float a[4], const float x[2])
{
//!# b=a'*x
	b[0] = 0.0; b[1] = 0.0;
	mc_addmulatx2x2ff(b, a, x);
}

MC_TARGET_FUNC void mc_mulatx2x2(double * b, const double a[4], const double x[2])
{
//!# b=a'*x
	const double x0 = x[0];
	const double x1 = x[1];

	b[0] = a[0] * x0;
	b[0] = b[0] + (a[2] * x1);

	b[1] = a[1] * x0;
	b[1] = b[1] + (a[3] * x1);
}

MC_TARGET_FUNC void mc_mulatx2x2l(long double * b, const long double a[4], const long double x[2])
{
//!# b=a'*x
	const long double x0 = x[0];
	const long double x1 = x[1];

	b[0] = a[0] * x0;
	b[0] = b[0] + (a[2] * x1);

	b[1] = a[1] * x0;
	b[1] = b[1] + (a[3] * x1);
}

#endif /* !MC_MULATX2X2_H */

/* EOF */