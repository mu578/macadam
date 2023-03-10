//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_intge.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_twosum.h>

#ifndef MC_INTGE_H
#define MC_INTGE_H

#pragma mark - mc_intge -

MC_TARGET_PROC float mc_intgef(float (*fx)(const float x), const float a, const float b, const int n)
{
//!# Computing the integral of f(x) within the interval [a, b].
	int i   = 1;
	float e = 0.0f, s = 0.0f, y;
	if (n > 0) {
		if (n < 0x1000001) {
			const float dt = (b - a) / mc_cast(const float, n);
			s              = fx(a + 0.5f * dt);
			for (; i < n;  i++) {
				mc_twosumf(s, fx(a + (i + 0.5f) * dt), &s, &y);
				e  = e + y;
			}
			return (s + e) * dt;
		} else {
			const double dt = mc_cast_expr(const double, (b - a)) / mc_cast(const double, n);
			s               = fx(a + 0.5f * mc_cast(const float, dt));
			for (; i < n;  i++) {
				mc_twosumf(s, fx(a + (i + 0.5f) * mc_cast(const float, dt)), &s, &y);
				e  = e + y;
			}
			return (s + e) * mc_cast(const float, dt);
		}
	}
	return MCK_NAN;
}

MC_TARGET_PROC double mc_intgeff(float (*fx)(const float x), const float a, const float b, const int n)
{
//!# Computing the integral of f(x) within the interval [a, b].
	const double ad   = mc_cast(double, a);
	const double bd   = mc_cast(double, b);
	int i             = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		const double dt = (bd - ad) / mc_cast(const double, n);
		s               = mc_cast_expr(double, fx(mc_cast_expr(float, ad + 0.5 * dt)));
		for (; i < n; i++) {
			mc_twosum(s, mc_cast_expr(double, fx(mc_cast_expr(float, ad + (i + 0.5) * dt))), &s, &y);
			e  = e + y;
		}
		return (s + e) * dt;
	}
	return MCK_NAN;
}

MC_TARGET_PROC double mc_intge(double (*fx)(const double x), const double a, const double b, const int n)
{
//!# Computing the integral of f(x) within the interval [a, b].
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		const double dt = (b - a) / mc_cast(const double, n);
		s               = fx(a + 0.5 * dt);
		for (; i < n; i++) {
			mc_twosum(s, fx(a + (i + 0.5) * dt), &s, &y);
			e  = e + y;
		}
		return (s + e) * dt;
	}
	return MCK_NAN;
}

MC_TARGET_PROC long double mc_intgel(long double (*fx)(const long double x), const long double a, const long double b, const int n)
{
//!# Computing the integral of f(x) within the interval [a, b].
	int i         = 1;
	long double e = 0.0L, s = 0.0L, y;
	if (n > 0) {
		const long double dt = (b - a) / mc_cast(const long double, n);
		s                    = fx(a + 0.5L * dt);
		for (; i < n; i++) {
			mc_twosuml(s, fx(a + (i + 0.5L) * dt), &s, &y);
			e  = e + y;
		}
		return (s + e) * dt;
	}
	return MCK_NAN;
}

#endif /* !MC_INTGE_H */

/* EOF */