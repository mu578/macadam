//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_outpxy3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_OUTPXY3X3_H
#define MC_OUTPXY3X3_H

#pragma mark - mc_outpxy3x3 -

MC_TARGET_FUNC void mc_outpxy3x3f(float a[9], const float x[3], const float y[3])
{
//!# Requires a[3 x 3], x[3 x 1] and y[3 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * y[0];
	a[1] = x[0] * y[1];
	a[2] = x[0] * y[2];

	a[3] = x[1] * y[0];
	a[4] = x[1] * y[1];
	a[5] = x[1] * y[2];

	a[6] = x[2] * y[0];
	a[7] = x[2] * y[1];
	a[8] = x[2] * y[2];
}

MC_TARGET_FUNC void mc_outpxy3x3ff(double a[9], const float x[3], const float y[3])
{
//!# Requires a[3 x 3], x[3 x 1] and y[3 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = mc_cast(double, x[0]) * mc_cast(double, y[0]);
	a[1] = mc_cast(double, x[0]) * mc_cast(double, y[1]);
	a[2] = mc_cast(double, x[0]) * mc_cast(double, y[2]);

	a[3] = mc_cast(double, x[1]) * mc_cast(double, y[0]);
	a[4] = mc_cast(double, x[1]) * mc_cast(double, y[1]);
	a[5] = mc_cast(double, x[1]) * mc_cast(double, y[2]);

	a[6] = mc_cast(double, x[2]) * mc_cast(double, y[0]);
	a[7] = mc_cast(double, x[2]) * mc_cast(double, y[1]);
	a[8] = mc_cast(double, x[2]) * mc_cast(double, y[2]);
}

MC_TARGET_FUNC void mc_outpxy3x3fd(double a[9], const float x[3], const double y[3])
{
//!# Requires a[3 x 3], x[3 x 1] and y[3 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = mc_cast(double, x[0]) * y[0];
	a[1] = mc_cast(double, x[0]) * y[1];
	a[2] = mc_cast(double, x[0]) * y[2];

	a[3] = mc_cast(double, x[1]) * y[0];
	a[4] = mc_cast(double, x[1]) * y[1];
	a[5] = mc_cast(double, x[1]) * y[2];

	a[6] = mc_cast(double, x[2]) * y[0];
	a[7] = mc_cast(double, x[2]) * y[1];
	a[8] = mc_cast(double, x[2]) * y[2];
}

MC_TARGET_FUNC void mc_outpxy3x3df(double a[9], const double x[3], const float y[3])
{
//!# Requires a[3 x 3], x[3 x 1] and y[3 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * mc_cast(double, y[0]);
	a[1] = x[0] * mc_cast(double, y[1]);
	a[2] = x[0] * mc_cast(double, y[2]);

	a[3] = x[1] * mc_cast(double, y[0]);
	a[4] = x[1] * mc_cast(double, y[1]);
	a[5] = x[1] * mc_cast(double, y[2]);

	a[6] = x[2] * mc_cast(double, y[0]);
	a[7] = x[2] * mc_cast(double, y[1]);
	a[8] = x[2] * mc_cast(double, y[2]);
}

MC_TARGET_FUNC void mc_outpxy3x3(double a[9], const double x[3], const double y[3])
{
//!# Requires a[3 x 3], x[3 x 1] and y[3 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * y[0];
	a[1] = x[0] * y[1];
	a[2] = x[0] * y[2];

	a[3] = x[1] * y[0];
	a[4] = x[1] * y[1];
	a[5] = x[1] * y[2];

	a[6] = x[2] * y[0];
	a[7] = x[2] * y[1];
	a[8] = x[2] * y[2];
}

MC_TARGET_FUNC void mc_outpxy3x3l(long double a[9], const long double x[3], const long double y[3])
{
//!# Requires a[3 x 3], x[3 x 1] and y[3 x 1].
//!# c=x*y' i.e outer product of two vectors.
	a[0] = x[0] * y[0];
	a[1] = x[0] * y[1];
	a[2] = x[0] * y[2];

	a[3] = x[1] * y[0];
	a[4] = x[1] * y[1];
	a[5] = x[1] * y[2];

	a[6] = x[2] * y[0];
	a[7] = x[2] * y[1];
	a[8] = x[2] * y[2];
}

#endif /* !MC_OUTPXY3X3_H */

/* EOF */