//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_hurwitz_zeta.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_exp2m1.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_lchoose.h>
#include <macadam/details/math/mc_pow.h>
#include <macadam/details/math/mc_riemann_zeta.h>

#ifndef MC_HURWITZ_ZETA_H
#define MC_HURWITZ_ZETA_H

#pragma mark - mc_hurwitz_zeta -

MC_TARGET_FUNC float mc_hurwitz_zetaf(const float s, const float q)
{
	const float w  = MCK_KF(FLT_MAX_10_EXP) * MCK_KF(MCK_LOGE10);
	float r        = 0.0f, c, term, sign;
	unsigned int i = 0, j;

	if ((mc_isnan(s) || mc_isinf(s)) || (mc_isnan(q) || mc_isinf(q)) || (s == 1.0f || q <= 0.0f)) {
		return MCK_NAN;
	}
	if (s == 0.0f) {
		return 0.5f - q;
	}
	if (q == 0.5f && s < 16384.0f) {
		return mc_exp2m1f(s) * mc_riemann_zetaf(s);
	}
	if (q == 1.0f) {
		return mc_riemann_zetaf(s);
	}
	for (; i < 10000; i++) {
		sign = 1.0f;
		term = 0.0f;
		for (j = 0; j <= i; j++) {
			if (w < (c = mc_lchoosef(mc_cast(const float, i), mc_cast(const float, j)))) {
				c    = sign * mc_expf(c);
				sign = sign * -1.0f;
				term = term + (c * mc_powf(q + mc_cast(const float, j), -s));
				continue;
			}
			goto hurwitz_end;
		}
		term = term * (1.0f / mc_cast(const float, (i + 1)));
		if (mc_fabsf(term * (1.0f / r)) < MCLIMITS_EPSILONF) {
			goto hurwitz_end;
		}
		r = r + term;
	}
hurwitz_end:
	return r / (s - 1.0f);
}

MC_TARGET_FUNC double mc_hurwitz_zeta(const double s, const double q)
{
	const double w = MCK_K(DBL_MAX_10_EXP) * MCK_K(MCK_LOGE10);
	double r       = 0.0, c, term, sign;
	unsigned int i = 0, j;

	if ((mc_isnan(s) || mc_isinf(s)) || (mc_isnan(q) || mc_isinf(q)) || (s == 1.0 || q <= 0.0)) {
		return MCK_NAN;
	}
	if (s == 0.0) {
		return 0.5 - q;
	}
	if (q == 0.5 && s < 16384.0) {
		return mc_exp2m1(s) * mc_riemann_zeta(s);
	}
	if (q == 1.0) {
		return mc_riemann_zeta(s);
	}
	for (; i < 10000; i++) {
		sign = 1.0;
		term = 0.0;
		for (j = 0; j <= i; j++) {
			if (w < (c = mc_lchoose(mc_cast(const double, i), mc_cast(const double, j)))) {
				c    = sign * mc_exp(c);
				sign = sign * -1.0;
				term = term + (c * mc_pow(q + mc_cast(const double, j), -s));
				continue;
			}
			goto hurwitz_end;
		}
		term = term * (1.0 / mc_cast(const double, (i + 1)));
		if (mc_fabs(term * (1.0 / r)) < MCLIMITS_EPSILON) {
			goto hurwitz_end;
		}
		r = r + term;
	}
hurwitz_end:
	return r / (s - 1.0);
}

MC_TARGET_FUNC long double mc_hurwitz_zetal(const long double s, const long double q)
{
	const long double w = MCK_KL(LDBL_MAX_10_EXP) * MCK_KL(MCK_LOGE10);
	long double r       = 0.0L, c, term, sign;
	unsigned int i       = 0, j;

	if ((mc_isnan(s) || mc_isinf(s)) || (mc_isnan(q) || mc_isinf(q)) || (s == 1.0L || q <= 0.0L)) {
		return MCK_NAN;
	}
	if (s == 0.0L) {
		return 0.5L - q;
	}
	if (q == 0.5L && s < 16384.0L) {
		return mc_exp2m1l(s) * mc_riemann_zetal(s);
	}
	if (q == 1.0L) {
		return mc_riemann_zetal(s);
	}
	for (; i < 10000; i++) {
		sign = 1.0L;
		term = 0.0L;
		for (j = 0; j <= i; j++) {
			if (w < (c = mc_lchoosel(mc_cast(const long double, i), mc_cast(const long double, j)))) {
				c    = sign * mc_expl(c);
				sign = sign * -1.0L;
				term = term + (c * mc_powl(q + mc_cast(const long double, j), -s));
				continue;
			}
			goto hurwitz_end;
		}
		term = term * (1.0L / mc_cast(const long double, (i + 1)));
		if (mc_fabsl(term * (1.0L / r)) < MCLIMITS_EPSILONL) {
			goto hurwitz_end;
		}
		r = r + term;
	}
hurwitz_end:
	return r / (s - 1.0L);
}

#endif /* !MC_HURWITZ_ZETA_H */

/* EOF */