//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_legendre_pn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_raise3.h>
#include <macadam/details/math/mc_raise4.h>
#include <macadam/details/math/mc_raise5.h>

#ifndef MC_LEGENDRE_PN_H
#define MC_LEGENDRE_PN_H

#pragma mark - mc_legendre_p0 -

MC_TARGET_PROC float mc_legendre_p0f(const float x)
{
	mc_unused(x);
	return 1.0f;
}

MC_TARGET_PROC double mc_legendre_p0(const double x)
{
	mc_unused(x);
	return 1.0;
}

MC_TARGET_PROC long double mc_legendre_p0l(const long double x)
{
	mc_unused(x);
	return 1.0L;
}

#pragma mark - mc_legendre_p1 -

MC_TARGET_PROC float mc_legendre_p1f(const float x)
{
	return x;
}

MC_TARGET_PROC double mc_legendre_p1(const double x)
{
	return x;
}

MC_TARGET_PROC long double mc_legendre_p1l(const long double x)
{
	return x;
}

#pragma mark - mc_legendre_p2 -

MC_TARGET_PROC float mc_legendre_p2f(const float x)
{
	return 1.5f * mc_raise2f(x) - 0.5f;
}

MC_TARGET_PROC double mc_legendre_p2(const double x)
{
	return 1.5 * mc_raise2(x) - 0.5;
}

MC_TARGET_PROC long double mc_legendre_p2l(const long double x)
{
	return 1.5L * mc_raise2l(x) - 0.5L;
}

#pragma mark - mc_legendre_p3 -

MC_TARGET_PROC float mc_legendre_p3f(const float x)
{
	return 0.5f * ((5.0f * mc_raise3f(x)) - (3.0f * x));
}

MC_TARGET_PROC double mc_legendre_p3(const double x)
{
	return 0.5 * ((5.0 * mc_raise3(x)) - (3.0 * x));
}

MC_TARGET_PROC long double mc_legendre_p3l(const long double x)
{
	return 0.5L * ((5.0L * mc_raise3l(x)) - (3.0L * x));
}

#pragma mark - mc_legendre_p4 -

MC_TARGET_PROC float mc_legendre_p4f(const float x)
{
	return 0.125f * (((35.0f * mc_raise4f(x)) - (30.0f * mc_raise2f(x))) + 3.0f);
}

MC_TARGET_PROC double mc_legendre_p4(const double x)
{
	return 0.125 * (((35.0 * mc_raise4(x)) - (30.0 * mc_raise2(x))) + 3.0);
}

MC_TARGET_PROC long double mc_legendre_p4l(const long double x)
{
	return 0.125L * (((35.0L * mc_raise4l(x)) - (30.0L * mc_raise2l(x))) + 3.0L);
}

#pragma mark - mc_legendre_p5 -

MC_TARGET_PROC float mc_legendre_p5f(const float x)
{
	return 0.125f * (((63.0f * mc_raise5f(x)) - (70.0f * mc_raise3f(x))) + (15.0f * x));
}

MC_TARGET_PROC double mc_legendre_p5(const double x)
{
	return 0.125 * (((63.0 * mc_raise5(x)) - (70.0 * mc_raise3(x))) + (15.0 * x));
}

MC_TARGET_PROC long double mc_legendre_p5l(const long double x)
{
	return 0.125L * (((63.0L * mc_raise5l(x)) - (70.0L * mc_raise3l(x))) + (15.0L * x));
}

#pragma mark - mc_legendre_p6 -

MC_TARGET_PROC float mc_legendre_p6f(const float x)
{
	return 0.0625f * (((231.0f * mc_raise6f(x)) - (315.0f * mc_raise4f(x)) + (105.0f * mc_raise2f(x))) - 5.0f);
}

MC_TARGET_PROC double mc_legendre_p6(const double x)
{
	return 0.0625 * (((231.0 * mc_raise6(x)) - (315.0 * mc_raise4(x)) + (105.0 * mc_raise2(x))) - 5.0);
}

MC_TARGET_PROC long double mc_legendre_p6l(const long double x)
{
	return 0.0625L * (((231.0L * mc_raise6l(x)) - (315.0L * mc_raise4l(x)) + (105.0L * mc_raise2l(x))) - 5.0L);
}

#pragma mark - mc_legendre_pn -

MC_TARGET_FUNC float mc_legendre_pnf(const unsigned int n, const float x)
{
//!# Legendre polynomials or functions.
	float pli      = 0.0f;
	float pl0      = 0.0f;
	float pl1      = 0.0f;
	unsigned int i = 2;
	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	if (x == 1.0f || x == -1.0f) {
		pl1 = 1.0f;
		if ((x < 0.0f) && ((n % 2) != 0)) {
			pl1 = -(pl1);
		}
	} else if (n >= 1 && n < MCLIMITS_IMAX) {
		if (n == 1) {
			pl1 = x;
		} else if (n == 2) {
			pl1 = mc_legendre_p2f(x);
		} else if ((x == 0.0f) && ((n % 2) != 0)) {
			pl1 = 0.0f;
		} else {
			pl0 = 1.0f;
			pl1 = x;
			for (; i <= n; i++) {
				pli = ((2.0f * mc_cast(float, i) - 1.0f) * x * pl1 - mc_cast(float, (i - 1)) * pl0) / mc_cast(const float, i);
				pl0 = pl1;
				pl1 = pli;
			}
		}
	}
	return pl1;
}

MC_TARGET_FUNC double mc_legendre_pn(const unsigned int n, const double x)
{
//!# Legendre polynomials or functions.
	double pli     = 0.0;
	double pl0     = 0.0;
	double pl1     = 0.0;
	unsigned int i = 2;
	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	if (x == 1.0 || x == -1.0) {
		pl1 = 1.0;
		if ((x < 0.0) && ((n % 2) != 0)) {
			pl1 = -(pl1);
		}
	} else if (n >= 1 && n < MCLIMITS_IMAX) {
		if (n == 1) {
			pl1 = x;
		} else if (n == 2) {
			pl1 = mc_legendre_p2(x);
		} else if ((x == 0.0) && ((n % 2) != 0)) {
			pl1 = 0.0;
		} else {
			pl0 = 1.0;
			pl1 = x;
			for (; i <= n; i++) {
				pli = ((2.0 * mc_cast(double, i) - 1.0) * x * pl1 - mc_cast(double, (i - 1)) * pl0) / mc_cast(const double, i);
				pl0 = pl1;
				pl1 = pli;
			}
		}
	}
	return pl1;
}

MC_TARGET_FUNC long double mc_legendre_pnl(const unsigned int n, const long double x)
{
//!# Legendre polynomials or functions.
	long double pli = 0.0L;
	long double pl0 = 0.0L;
	long double pl1 = 0.0L;
	unsigned int i  = 2;
	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	if (x == 1.0L || x == -1.0L) {
		pl1 = 1.0L;
		if ((x < 0.0L) && ((n % 2) != 0)) {
			pl1 = -(pl1);
		}
	} else if (n >= 1 && n < MCLIMITS_IMAX) {
		if (n == 1) {
			pl1 = x;
		} else if (n == 2) {
			pl1 = mc_legendre_p2l(x);
		} else if ((x == 0.0L) && ((n % 2) != 0)) {
			pl1 = 0.0L;
		} else {
			pl0 = 1.0L;
			pl1 = x;
			for (; i <= n; i++) {
				pli = ((2.0L * mc_cast(long double, i) - 1.0L) * x * pl1 - mc_cast(long double, (i - 1)) * pl0) / mc_cast(const long double, i);
				pl0 = pl1;
				pl1 = pli;
			}
		}
	}
	return pl1;
}

#endif /* !MC_LEGENDRE_PN_H */

/* EOF */