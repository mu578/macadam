//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_log1pmx.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_log1p.h>
#include <macadam/details/math/mc_logcf.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_LOG1PMX_H
#define MC_LOG1PMX_H

#pragma mark - mc_log1pmx -

MC_TARGET_FUNC float mc_log1pmxf(const float x)
{
//!# Computing log(1 + x) - x after Ian Smith's log1 function.
	const float tol = MCLIMITS_TINYF;

	float y , z, r = MCK_INFN;

	if (x != 1.0f) {
		if (x < -0.525f || x > 1.05f) {
			r = mc_log1pf(x) - x;
		} else {
			z = x / (2.0f + x);
			y = mc_raise2f(z);
			if (mc_fabsf(x) > 0.016f) {
				y = 2.0f * y * mc_logcff(y, 3.0f, 2.0f, tol);
			} else {
				y = ((((2.0f / 9.0f) * y + (2.0f / 7.0f)) * y + (2.0f / 5.0f)) * y + (2.0f / 3.0f)) * y;
			}
			r = z * (y - x);
		}
	}
	return r;
}

MC_TARGET_FUNC double mc_log1pmx(const double x)
{
//!# Computing log(1 + x) - x after Ian Smith's log1 function.
	const double tol = MCLIMITS_TINY;

	double y , z, r = MCK_INFN;
	if (x != 1.0) {
		if (x < -0.525 || x > 1.05) {
			r = mc_log1p(x) - x;
		} else {
			z = x / (2.0 + x);
			y = mc_raise2(z);
			if (mc_fabs(x) > 0.016) {
				y = 2.0 * y * mc_logcf(y, 3.0, 2.0, tol);
			} else {
				y = ((((2.0 / 9.0) * y + (2.0 / 7.0)) * y + (2.0 / 5.0)) * y + (2.0 / 3.0)) * y;
			}
			r = z * (y - x);
		}
	}
	return r;
}

MC_TARGET_FUNC long double mc_log1pmxl(const long double x)
{
//!# Computing log(1 + x) - x after Ian Smith's log1 function.
#	if MC_TARGET_HAVE_LONG_DOUBLE
	const long double tol = MCLIMITS_TINYL;
#	else
	const long double tol = MCK_KL(MCLIMITS_TINY);
#	endif

	long double y , z, r = MCK_INFN;
	if (x != 1.0L) {
		if (x < -0.525L || x > 1.05L) {
			r = mc_log1pl(x) - x;
		} else {
			z = x / (2.0L + x);
			y = mc_raise2l(z);
			if (mc_fabsl(x) > 0.016L) {
				y = 2.0L * y * mc_logcfl(y, 3.0L, 2.0L, tol);
			} else {
				y = ((((2.0L / 9.0L) * y + (2.0L / 7.0L)) * y + (2.0L / 5.0L)) * y + (2.0L / 3.0L)) * y;
			}
			r = z * (y - x);
		}
	}
	return r;
}

#endif /* !MC_LOG1PMX_H */

/* EOF */