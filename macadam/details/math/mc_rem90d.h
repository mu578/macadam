//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rem90d.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fisval.h>
#include <macadam/details/math/mc_floor.h>
#include <macadam/details/math/mc_fmod.h>
#include <macadam/details/math/mc_itrunc64.h>

#ifndef MC_REM90D_H
#define MC_REM90D_H

#pragma mark - mc_rem90d -

MC_TARGET_PROC int64_t mc_rem90df(const float x, float * z)
{
//!# z = x mod 90 argument reduction similar to @remint2 for degree given
//!# angles. Reducing x in degrees mod 90; z = x mod 90, |z| < 45. @todo
//!# improving situation in @sind and @cosd respectively.
	int64_t r  = 0;
	float w, y = 0.0f;

	if (mc_fisvalf(x) && x != 0.0f) {
		w = x;
		//!# argument is out of any computable range, accurency
		//!# is already lost, reducing first over a 360 field.
		if (mc_fabsf(x) >= 5826.0f) {
			w = mc_fmodf(w, 360.0f);
		}
		y = mc_floorf(w / 45.0f);
		r = mc_itrunc64f(y - 16.0f * mc_floorf(y / 16.0f));
		if (r & 1) {
			++r;
			y = y + 1.0f;
		}
		r = (mc_cast(uint64_t, r) >> 1) & 7;
		y = w - y * 45.0f;
	}
	*z = y;
	return r;
}

MC_TARGET_PROC int64_t mc_rem90d(const double x, double * z)
{
//!# z = x mod 90 argument reduction similar to @remint2 for degree given
//!# angles. Reducing x in degrees mod 90; z = x mod 90, |z| < 45. @todo
//!# improving situation in @sind and @cosd respectively.
	int64_t r   = 0;
	double w, y = 0.0;

	if (mc_fisval(x) && x != 0.0) {
		w = x;
		//!# argument is out of any computable range, accurency
		//!# is already lost, reducing first over a 360 field.
		if (mc_fabs(x) >= 23860928.0) {
			w = mc_fmod(w, 360.0);
		}
		y = mc_floor(w / 45.0);
		r = mc_itrunc64(y - 16.0 * mc_floor(y / 16.0));
		if (r & 1) {
			++r;
			y = y + 1.0;
		}
		r = (mc_cast(uint64_t, r) >> 1) & 7;
		y = w - y * 45.0;
	}
	*z = y;
	return r;
}

MC_TARGET_PROC int64_t mc_rem90dl(const long double x, long double * z)
{
//!# z = x mod 90 argument reduction similar to @remint2 for degree given
//!# angles. Reducing x in degrees mod 90; z = x mod 90, |z| < 45. @todo
//!# improving situation in @sind and @cosd respectively.
	int64_t r        = 0;
	long double w, y = 0.0L;

	if (mc_fisvall(x) && x != 0.0L) {
		w = x;
		//!# argument is out of any computable range, accurency
		//!# is already lost, reducing first over a 360 field.
		if (mc_fabsl(x) >= 3054198966.0L) {
			w = mc_fmodl(w, 360.0L);
		}
		y = mc_floorl(w / 45.0L);
		r = mc_itrunc64l(y - 16.0L * mc_floorl(y / 16.0L));
		if (r & 1) {
			++r;
			y = y + 1.0L;
		}
		r = (mc_cast(uint64_t, r) >> 1) & 7;
		y = w - y * 45.0L;
	}
	*z = y;
	return r;
}

#endif /* !MC_REM90D_H */

/* EOF */