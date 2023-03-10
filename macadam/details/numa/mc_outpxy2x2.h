//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_outpxy2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_OUTPXY2X2_H
#define MC_OUTPXY2X2_H

#pragma mark - mc_outpxy2x2 -

MC_TARGET_FUNC void mc_outpxy2x2f(float a[4], const float x[2], const float y[2])
{
//!# Requires a[2 x 2], x[2 x 1] and y[2 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * y[0];
	a[1] = x[0] * y[1];

	a[2] = x[1] * y[0];
	a[2] = x[1] * y[1];
}

MC_TARGET_FUNC void mc_outpxy2x2ff(double a[4], const float x[2], const float y[2])
{
//!# Requires a[2 x 2], x[2 x 1] and y[2 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = mc_cast(double, x[0]) * mc_cast(double, y[0]);
	a[1] = mc_cast(double, x[0]) * mc_cast(double, y[1]);

	a[2] = mc_cast(double, x[1]) * mc_cast(double, y[0]);
	a[2] = mc_cast(double, x[1]) * mc_cast(double, y[1]);
}

MC_TARGET_FUNC void mc_outpxy2x2fd(double a[4], const float x[2], const double y[2])
{
//!# Requires a[2 x 2], x[2 x 1] and y[2 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = mc_cast(double, x[0]) * y[0];
	a[1] = mc_cast(double, x[0]) * y[1];

	a[2] = mc_cast(double, x[1]) * y[0];
	a[2] = mc_cast(double, x[1]) * y[1];
}

MC_TARGET_FUNC void mc_outpxy2x2df(double a[4], const double x[2], const float y[2])
{
//!# Requires a[2 x 2], x[2 x 1] and y[2 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * mc_cast(double, y[0]);
	a[1] = x[0] * mc_cast(double, y[1]);

	a[2] = x[1] * mc_cast(double, y[0]);
	a[2] = x[1] * mc_cast(double, y[1]);
}

MC_TARGET_FUNC void mc_outpxy2x2(double a[4], const double x[2], const double y[2])
{
//!# Requires a[2 x 2], x[2 x 1] and y[2 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * y[0];
	a[1] = x[0] * y[1];

	a[2] = x[1] * y[0];
	a[2] = x[1] * y[1];
}

MC_TARGET_FUNC void mc_outpxy2x2l(long double a[4], const long double x[2], const long double y[2])
{
//!# Requires a[2 x 2], x[2 x 1] and y[2 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * y[0];
	a[1] = x[0] * y[1];

	a[2] = x[1] * y[0];
	a[2] = x[1] * y[1];
}

#endif /* !MC_OUTPXY2X2_H */

/* EOF */