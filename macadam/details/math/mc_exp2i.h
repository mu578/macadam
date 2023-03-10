//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_exp2i.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_absmag.h>
#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_isnormal.h>

#ifndef MC_EXP2I_H
#define MC_EXP2I_H

#pragma mark - mc_exp2i -

MC_TARGET_FUNC float mc_exp2if(const int e)
{
	float r       = 1.0f;
	int i         = 0;
	const float s = e >= 0 ? 2.0f : 0.5f;
	for (; i < mc_iabs(e); i++) {
		r = r * s;
		if (!mc_isfinite(r) || !mc_isnormal(r)) {
			r = e >= 0 ? MCK_INF : 0.0f;
			break;
		}
	}
	return r;
}

MC_TARGET_FUNC double mc_exp2i(const int e)
{
	double r       = 1.0;
	int i          = 0;
	const double s = e >= 0 ? 2.0 : 0.5;
	for (; i < mc_iabs(e); i++) {
		r = r * s;
		if (!mc_isfinite(r) || !mc_isnormal(r)) {
			r = e >= 0 ? MCK_INF : 0.0;
			break;
		}
	}
	return r;
}

MC_TARGET_FUNC long double mc_exp2il(const int e)
{
	long double r       = 1.0L;
	int i               = 0;
	const long double s = e >= 0 ? 2.0L : 0.5L;
	for (; i < mc_iabs(e); i++) {
		r = r * s;
		if (!mc_isfinite(r) || !mc_isnormal(r)) {
			r = e >= 0 ? MCK_INF : 0.0L;
			break;
		}
	}
	return r;
}

#endif /* !MC_EXP2I_H */

/* EOF */