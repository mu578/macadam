//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_exp10i.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_absmag.h>
#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_isnormal.h>

#ifndef MC_EXP10I_H
#define MC_EXP10I_H

#pragma mark - mc_exp10i -

MC_TARGET_FUNC float mc_exp10if(const int e)
{
	long double r       = 1.0L;
	float x             = 1.0f;
	int i               = 0;
	const long double s = e >= 0 ? 10.0L : 0.1L;
	for (; i < mc_iabs(e); i++) {
		r = (r * s);
		x = mc_cast(float, r);
		if (!mc_isfinite(x) || !mc_isnormal(x)) {
			x = e >= 0 ? MCK_INF : 0.0f;
			break;
		}
	}
	return x;
}

MC_TARGET_FUNC double mc_exp10i(const int e)
{
	long double r       = 1.0L;
	double x            = 1.0;
	int i               = 0;
	const long double s = e >= 0 ? 10.0L : 0.1L;
	for (; i < mc_iabs(e); i++) {
		r = (r * s);
		x = mc_cast(double, r);
		if (!mc_isfinite(x) || !mc_isnormal(x)) {
			x = e >= 0 ? MCK_INF : 0.0;
			break;
		}
	}
	return x;
}

MC_TARGET_FUNC long double mc_exp10il(const int e)
{
	long double r       = 1.0L;
	int i               = 0;
	const long double s = e >= 0 ? 10.0L : 0.1L;
	for (; i < mc_iabs(e); i++) {
		r = (r * s);
		if (!mc_isfinite(r) || !mc_isnormal(r)) {
			r = e >= 0 ? MCK_INF : 0.0L;
			break;
		}
	}
	return r;
}

#endif /* !MC_EXP10I_H */

/* EOF */