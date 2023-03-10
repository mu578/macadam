//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lmgamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gammaln.h>

#ifndef MC_LMGAMMA_H
#define MC_LMGAMMA_H

MC_TARGET_FUNC float mc_lmgammaf(unsigned int p, float a)
{
//!# Log of the multivariate gamma function.
	if (p == 0) {
		return MCK_INFP;
	}
	if (p == 1) {
		return mc_gammalnf(a);
	} else if (p < MCLIMITS_IMAX) {
		const float d  = mc_cast(float, p);
		float x        = d * (d - 1.0f) / 4.0f * MCK_KF(MCK_LOGEPI);
		unsigned int k = 0;
		for (; k < p; k++) {
			x = x + mc_gammalnf(a - 0.5f * mc_cast(float, k));
		}
		return x;
	}
	return MCK_NAN;
}

MC_TARGET_FUNC double mc_lmgamma(unsigned int p, double a)
{
//!# Log of the multivariate gamma function.
	if (p == 0) {
		return MCK_INFP;
	}
	if (p == 1) {
		return mc_gammaln(a);
	} else if (p < MCLIMITS_IMAX) {
		const double d = mc_cast(double, p);
		double x       = d * (d - 1.0) / 4.0 * MCK_K(MCK_LOGEPI);
		unsigned int k = 0;
		for (; k < p; k++) {
			x = x + mc_gammaln(a - 0.5 * mc_cast(double, k));
		}
		return x;
	}
	return MCK_NAN;
}

MC_TARGET_FUNC long double mc_lmgammal(unsigned int p, long double a)
{
//!# Log of the multivariate gamma function.
	if (p == 0) {
		return MCK_INFP;
	}
	if (p == 1) {
		return mc_gammalnl(a);
	} else if (p < MCLIMITS_IMAX) {
		const long double d = mc_cast(long double, p);
		long double x       = d * (d - 1.0L) / 4.0L * MCK_KL(MCK_LOGEPI);
		unsigned int k      = 0;
		for (; k < p; k++) {
			x = x + mc_gammalnl(a - 0.5L * mc_cast(long double, k));
		}
		return x;
	}
	return MCK_NAN;
}

#endif /* !MC_LMGAMMA_H */

/* EOF */