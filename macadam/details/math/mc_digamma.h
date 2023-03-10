//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_digamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cot.h>
#include <macadam/details/math/mc_evalpoly.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_pow.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_DIGAMMA_H
#define MC_DIGAMMA_H

#pragma mark - mc_digamma_approx0 -

MC_TARGET_PROC float mc_digammaf_approx0(const float x)
{
//!# Returns psi0(x) by generalized power series expansion.
	if (x <= 1E-5) {
		return -MCK_KF(MCK_G) - (1.0f / x);
	}
	if (x < 8.5f) {
		return mc_digammaf_approx0(1.0f + x) - 1.0f / x;
	}
	const float p = mc_powf(x, -2.0f);
	return (
		mc_logf(x) - 0.5f / x
		+ MCK_KF(MCK_ZETA_N1)  * p
		+ MCK_KF(MCK_ZETA_N3)  * mc_powf(p, 2)
		+ MCK_KF(MCK_ZETA_N5)  * mc_powf(p, 3)
		+ MCK_KF(MCK_ZETA_N7)  * mc_powf(p, 4)
		+ MCK_KF(MCK_ZETA_N9)  * mc_powf(p, 5)
		+ MCK_KF(MCK_ZETA_N11) * mc_powf(p, 6)
	);
}

MC_TARGET_PROC double mc_digamma_approx0(const double x)
{
//!# Returns psi0(x) by generalized power series expansion.
	if (x <= 1E-5) {
		return -MCK_K(MCK_G) - (1.0 / x);
	}
	if (x < 8.5) {
		return mc_digamma_approx0(1.0 + x) - 1.0 / x;
	}
	const double p = mc_pow(x, -2.0);
	return (
		mc_log(x) - 0.5 / x
		+ MCK_K(MCK_ZETA_N1)  * p
		+ MCK_K(MCK_ZETA_N3)  * mc_pow(p, 2)
		+ MCK_K(MCK_ZETA_N5)  * mc_pow(p, 3)
		+ MCK_K(MCK_ZETA_N7)  * mc_pow(p, 4)
		+ MCK_K(MCK_ZETA_N9)  * mc_pow(p, 5)
		+ MCK_K(MCK_ZETA_N11) * mc_pow(p, 6)
	);
}

MC_TARGET_PROC long double mc_digammal_approx0(const long double x)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
//!# Returns psi0(x) by generalized power series expansion.
	if (x <= 1E-5) {
		return -MCK_KL(MCK_G) - (1.0L / x);
	}
	if (x < 8.5L) {
		return mc_digammal_approx0(1.0L + x) - 1.0L / x;
	}
	const long double p = mc_powl(x, -2.0L);
	return (
		mc_logl(x) - 0.5L / x
		+ MCK_KL(MCK_ZETA_N1)  * p
		+ MCK_KL(MCK_ZETA_N3)  * mc_powl(p, 2)
		+ MCK_KL(MCK_ZETA_N5)  * mc_powl(p, 3)
		+ MCK_KL(MCK_ZETA_N7)  * mc_powl(p, 4)
		+ MCK_KL(MCK_ZETA_N9)  * mc_powl(p, 5)
		+ MCK_KL(MCK_ZETA_N11) * mc_powl(p, 6)
	);
#	else
	return mc_cast(long double, mc_digamma_approx0(mc_cast(const double, x)));
#	endif
}

#pragma mark - mc_digamma_approx1 -

MC_TARGET_PROC float mc_digammaf_approx1(const float x)
{
//!# Returns psi0(x) by Taylor series expansion.
	const float c1 = -8.33333333333333290000000000000000000000E-02f;
	const float c2 = +8.33333333333333320000000000000000000000E-03f;
	const float c3 = -3.96825396825396800000000000000000000000E-03f;
	const float c4 = +4.16666666666666660000000000000000000000E-03f;
	const float c5 = -7.57575757575757600000000000000000000000E-03f;
	const float c6 = +2.10927960927960940000000000000000000000E-02f;
	const float c7 = -8.33333333333333290000000000000000000000E-02f;

	float r = 0.0f, y, w = x;

	do {
		r = r - (1.0f / w);
		w = w + 1.0f;
	} while (w < 10.0f);

	y = mc_raise2f(w);
	r = r + (mc_logf(w) - 0.5f / w);
	w = y;

	r = r + (c1 * (1.0f / w));
	w = w * y;
	r = r + (c2 * (1.0f / w));
	w = w * y;
	r = r + (c3 * (1.0f / w));
	w = w * y;
	r = r + (c4 * (1.0f / w));
	w = w * y;
	r = r + (c5 * (1.0f / w));
	w = w * y;
	r = r + (c6 * (1.0f / w));
	w = w * y;
	r = r + (c7 * (1.0f / w));

	return r;
}

MC_TARGET_PROC double mc_digamma_approx1(const double x)
{
//!# Returns psi0(x) by Taylor series expansion.
	const double c1 = -8.3333333333333329000000000000000000000000E-02;
	const double c2 = +8.3333333333333332000000000000000000000000E-03;
	const double c3 = -3.9682539682539680000000000000000000000000E-03;
	const double c4 = +4.1666666666666666000000000000000000000000E-03;
	const double c5 = -7.5757575757575760000000000000000000000000E-03;
	const double c6 = +2.1092796092796094000000000000000000000000E-02;
	const double c7 = -8.3333333333333329000000000000000000000000E-02;

	double r = 0.0, y, w = x;

	do {
		r = r - (1.0 / w);
		w = w + 1.0;
	} while (w < 10.0);

	y = mc_raise2(w);
	r = r + (mc_log(w) - 0.5 / w);
	w = y;

	r = r + (c1 * (1.0 / w));
	w = w * y;
	r = r + (c2 * (1.0 / w));
	w = w * y;
	r = r + (c3 * (1.0 / w));
	w = w * y;
	r = r + (c4 * (1.0 / w));
	w = w * y;
	r = r + (c5 * (1.0 / w));
	w = w * y;
	r = r + (c6 * (1.0 / w));
	w = w * y;
	r = r + (c7 * (1.0 / w));

	return r;
}

MC_TARGET_PROC long double mc_digammal_approx1(const long double x)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
//!# Returns psi0(x) by Taylor series expansion.
	const long double c1 = -8.333333333333332900000000000000000000000000000000000000000000000E-02L;
	const long double c2 = +8.333333333333333200000000000000000000000000000000000000000000000E-03L;
	const long double c3 = -3.968253968253968000000000000000000000000000000000000000000000000E-03L;
	const long double c4 = +4.166666666666666600000000000000000000000000000000000000000000000E-03L;
	const long double c5 = -7.575757575757576000000000000000000000000000000000000000000000000E-03L;
	const long double c6 = +2.109279609279609400000000000000000000000000000000000000000000000E-02L;
	const long double c7 = -8.333333333333332900000000000000000000000000000000000000000000000E-02L;

	long double r = 0.0L, y, w = x;

	do {
		r = r - (1.0L / w);
		w = w + 1.0L;
	} while (w < 10.0L);

	y = mc_raise2l(w);
	r = r + (mc_logl(w) - 0.5L / w);
	w = y;

	r = r + (c1 * (1.0L / w));
	w = w * y;
	r = r + (c2 * (1.0L / w));
	w = w * y;
	r = r + (c3 * (1.0L / w));
	w = w * y;
	r = r + (c4 * (1.0L / w));
	w = w * y;
	r = r + (c5 * (1.0L / w));
	w = w * y;
	r = r + (c6 * (1.0L / w));
	w = w * y;
	r = r + (c7 * (1.0L / w));

	return r;
#	else
	return mc_cast(long double, mc_digamma_approx1(mc_cast(const double, x)));
#	endif
}

#pragma mark - mc_digamma_approx2 -

MC_TARGET_PROC float mc_digammaf_approx2(const float x)
{
//!# Lanczos expansion coefficients. @see zeta, b2n.
	const float C[] =
	{
		  +8.33333333333333290000000000000000000000E-02f
		, -8.33333333333333320000000000000000000000E-03f
		, +3.96825396825396800000000000000000000000E-03f
		, -4.16666666666666660000000000000000000000E-03f
		, +7.57575757575757600000000000000000000000E-03f
		, -2.10927960927960940000000000000000000000E-02f
		, +8.33333333333333290000000000000000000000E-02f
		, -4.43259803921568600000000000000000000000E-01f
	};

	float a = 0.0f, w, y, t;
	int i   = 1, n;

	w = y = x;
	if (w == 0.0f) {
		return MCK_INFP;
	}
	if (w < 0.0f) {
		a = -MCK_KF(MCK_PI) * mc_cotf(MCK_KF(MCK_PI) * y);
		y = 1 - y;
		w = y;
	}
	if (w < 7.0f) {
		n =  7 - mc_cast(int, mc_floorf(w));
		for (; i < n; i++) {
			a = a - (1.0f / (y + mc_cast(float, i)));
		}
		a = a - (1.0f / y);
		y = y + mc_cast(float, n);
	}
	t = 1.0f / y;
	a = mc_logf(y) - 0.5f * t;
	t = mc_raise2f(t);
	return a - (t * mc_evalpolyf(t, C, 8, 0));
}

MC_TARGET_PROC double mc_digamma_approx2(const double x)
{
//!# Lanczos expansion coefficients. @see zeta, b2n.
	const double C[] =
	{
		  +8.3333333333333329000000000000000000000000E-02
		, -8.3333333333333332000000000000000000000000E-03
		, +3.9682539682539680000000000000000000000000E-03
		, -4.1666666666666666000000000000000000000000E-03
		, +7.5757575757575760000000000000000000000000E-03
		, -2.1092796092796094000000000000000000000000E-02
		, +8.3333333333333329000000000000000000000000E-02
		, -4.4325980392156860000000000000000000000000E-01
	};

	double a = 0.0, w, y, t;
	int i    = 1, n;

	w = y = x;
	if (w == 0.0) {
		return MCK_INFP;
	}
	if (w < 0.0) {
		a = -MCK_K(MCK_PI) * mc_cot(MCK_K(MCK_PI) * y);
		y = 1 - y;
		w = y;
	}
	if (w < 7.0) {
		n =  7 - mc_cast(int, mc_floor(w));
		for (; i < n; i++) {
			a = a - (1.0 / (y + mc_cast(double, i)));
		}
		a = a - (1.0 / y);
		y = y + mc_cast(double, n);
	}
	t = 1.0 / y;
	a = mc_log(y) - 0.5 * t;
	t = mc_raise2(t);
	return a - (t * mc_evalpoly(t, C, 8, 0));
}

MC_TARGET_PROC long double mc_digammal_approx2(const long double x)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
//!# Lanczos expansion coefficients. @see zeta, b2n.
	const long double C[] =
	{
		  +8.333333333333332900000000000000000000000000000000000000000000000E-02L
		, -8.333333333333333200000000000000000000000000000000000000000000000E-03L
		, +3.968253968253968000000000000000000000000000000000000000000000000E-03L
		, -4.166666666666666600000000000000000000000000000000000000000000000E-03L
		, +7.575757575757576000000000000000000000000000000000000000000000000E-03L
		, -2.109279609279609400000000000000000000000000000000000000000000000E-02L
		, +8.333333333333332900000000000000000000000000000000000000000000000E-02L
		, -4.432598039215686000000000000000000000000000000000000000000000000E-01L
	};

	long double a = 0.0L, w, y, t;
	int i         = 1, n;

	w = y = x;
	if (w == 0.0L) {
		return MCK_INFP;
	}
	if (w < 0.0L) {
		a = -MCK_KL(MCK_PI) * mc_cotl(MCK_KL(MCK_PI) * y);
		y = 1 - y;
		w = y;
	}
	if (w < 7.0) {
		n =  7 - mc_cast(int, mc_floorl(w));
		for (; i < n; i++) {
			a = a - (1.0L / (y + mc_cast(long double, i)));
		}
		a = a - (1.0L / y);
		y = y + mc_cast(long double, n);
	}
	t = 1.0L / y;
	a = mc_logl(y) - 0.5L * t;
	t = mc_raise2l(t);
	return a - (t * mc_evalpolyl(t, C, 8, 0));
#	else
	return mc_cast(long double, mc_digamma_approx2(mc_cast(const double, x)));
#	endif
}

#pragma mark - mc_digamma -

MC_TARGET_FUNC float mc_digammaf(const float x)
{
	return mc_digammaf_approx2(x);
}

MC_TARGET_FUNC double mc_digamma(const double x)
{
	return mc_digamma_approx2(x);
}

MC_TARGET_FUNC long double mc_digammal(const long double x)
{
	return mc_digammal_approx2(x);
}

#endif /* !MC_DIGAMMA_H */

/* EOF */