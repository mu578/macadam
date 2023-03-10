//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_expi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_absmag.h>
#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_isnormal.h>

#ifndef MC_EXPI_H
#define MC_EXPI_H

#pragma mark - mc_expi -

MC_TARGET_FUNC float mc_expif(const float x, const int e)
{
	long double r = (e < 0
		? 1.0L / mc_cast(const long double, x)
		: mc_cast(const long double, x)
	);
	long double y = 1.0L;
	int         p = mc_iabs(e);
	if (p == 0) {
		return 1.0f;
	} else if (p == 1) {
		return mc_cast(float, r);
	} else {
		do {
			if (!(p & 1)) {
				r = mc_raise2l(r);
				p = p / 2;
			} else {
				y = r * y;
				r = mc_raise2l(r);
				p = (p - 1) / 2;
			}
			if (
				   !mc_isfinite(mc_cast_expr(float, r * y))
				|| !mc_isnormal(mc_cast_expr(float, r * y))
			) {
				r = e >= 0 ? MCK_INF : 0.0f;
				break;
			}
		} while (p > 1);
		return mc_cast_expr(float, r * y);
	}
	return MCK_NAN;
}

MC_TARGET_FUNC double mc_expi(const double x, const int e)
{
	long double r = (e < 0
		? 1.0L / mc_cast(const long double, x)
		: mc_cast(const long double, x)
	);
	long double y = 1.0L;
	int         p = mc_iabs(e);
	if (p == 0) {
		return 1.0;
	} else if (p == 1) {
		return mc_cast(double, r);
	} else {
		do {
			if (!(p & 1)) {
				r = mc_raise2l(r);
				p = p / 2;
			} else {
				y = r * y;
				r = mc_raise2l(r);
				p = (p - 1) / 2;
			}
			if (
				   !mc_isfinite(mc_cast_expr(double, r * y))
				|| !mc_isnormal(mc_cast_expr(double, r * y))
			) {
				r = e >= 0 ? MCK_INF : 0.0;
				break;
			}
		} while (p > 1);
		return mc_cast_expr(double, r * y);
	}
	return MCK_NAN;
}

MC_TARGET_FUNC long double mc_expil(const long double x, const int e)
{
	long double r = (e < 0
		? 1.0L / x
		: x
	);
	long double y = 1.0L;
	int         p = mc_iabs(e);
	if (p == 0) {
		return 1.0L;
	} else if (p == 1) {
		return r;
	} else {
		do {
			if (!(p & 1)) {
				r = mc_raise2l(r);
				p = p / 2;
			} else {
				y = r * y;
				r = mc_raise2l(r);
				p = (p - 1) / 2;
			}
			if (
				   !mc_isfinite((r * y))
				|| !mc_isnormal((r * y))
			) {
				r = e >= 0 ? MCK_INF : 0.0L;
				break;
			}
		} while (p > 1);
		return r * y;
	}
	return MCK_NAN;
}

#endif /* !MC_EXPI_H */

/* EOF */