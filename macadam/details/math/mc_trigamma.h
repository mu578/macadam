//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trigamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_floor.h>
#include <macadam/details/math/mc_rsqr.h>
#include <macadam/details/math/mc_sin.h>

#ifndef MC_TRIGAMMA_H
#define MC_TRIGAMMA_H

#pragma mark - mc_trigamma_approx0 -

MC_TARGET_PROC float mc_trigammaf_approx0(const float x)
{
	float g, y, z;
	if ((x <= 0) && (mc_floorf(x) == x)) {
		return MCK_INFP;
	}
	if ((x <= 0) && (mc_floorf(x) != x)) {
		return -mc_trigammaf_approx0(-x + 1.0f) + (MCK_KF(MCK_PI) / mc_sinf(-MCK_KF(MCK_PI) * x)) * (MCK_KF(MCK_PI) / mc_sinf(-MCK_KF(MCK_PI) * x));
	}
	if (x <= 1E-5f) {
		return mc_rsqrf(x);
	}
	z = x;
	g = 0.0f;
	while (z < 5.0f) {
		g = g + mc_rsqrf(z);
		z = z + 1.0f;
	}
	y  = mc_rsqrf(z);
//!# Expansion as a Laurent series.
	g = g + (0.5f * y + (1.0f + y * (MCK_KF(MCK_BN2) + y * (MCK_KF(MCK_BN4) + y * (MCK_KF(MCK_BN6) + y * MCK_KF(MCK_BN8))))) / z);
	return g;
}

MC_TARGET_PROC double mc_trigamma_approx0(const double x)
{
	double g, y, z;
	if ((x <= 0) && (mc_floor(x) == x)) {
		return MCK_INFP;
	}
	if ((x <= 0) && (mc_floor(x) != x)) {
		return -mc_trigamma_approx0(-x + 1.0) + (MCK_K(MCK_PI) / mc_sin(-MCK_K(MCK_PI) * x)) * (MCK_K(MCK_PI) / mc_sin(-MCK_K(MCK_PI) * x));
	}
	if (x <= 1E-5) {
		return mc_rsqr(x);
	}
	z = x;
	g = 0.0;
	while (z < 5.0) {
		g = g + mc_rsqr(z);
		z = z + 1.0;
	}
	y  = mc_rsqr(z);
//!# Expansion as a Laurent series.
	g = g + (0.5 * y + (1.0 + y * (MCK_K(MCK_BN2) + y * (MCK_K(MCK_BN4) + y * (MCK_K(MCK_BN6) + y * MCK_K(MCK_BN8))))) / z);
	return g;
}

MC_TARGET_PROC long double mc_trigammal_approx0(const long double x)
{
	long double g, y, z;
	if ((x <= 0) && (mc_floorl(x) == x)) {
		return MCK_INFP;
	}
	if ((x <= 0) && (mc_floorl(x) != x)) {
		return -mc_trigammal_approx0(-x + 1.0L) + (MCK_KL(MCK_PI) / mc_sinl(-MCK_KL(MCK_PI) * x)) * (MCK_KL(MCK_PI) / mc_sinl(-MCK_KL(MCK_PI) * x));
	}
	if (x <= 1E-5L) {
		return mc_rsqrl(x);
	}
	z = x;
	g = 0.0L;
	while (z < 5.0L) {
		g = g + mc_rsqrl(z);
		z = z + 1.0L;
	}
	y  = mc_rsqrl(z);
//!# Expansion as a Laurent series.
	g = g + (0.5L * y + (1.0L + y * (MCK_KL(MCK_BN2) + y * (MCK_KL(MCK_BN4) + y * (MCK_KL(MCK_BN6) + y * MCK_KL(MCK_BN8))))) / z);
	return g;
}

#pragma mark - mc_trigamma_approx1 -

MC_TARGET_PROC float mc_trigammaf_approx1(const float x)
{
//!# Returns psi1(x) by Taylor series expansion.
	const float c1 = +1.00000000000000000000000000000000000000E+00f;
	const float c2 = +1.66666666666666660000000000000000000000E-01f;
	const float c3 = -3.33333333333333330000000000000000000000E-02f;
	const float c4 = +2.38095238095238080000000000000000000000E-02f;
	const float c5 = -3.33333333333333330000000000000000000000E-02f;
	const float c6 = +7.57575757575757600000000000000000000000E-02f;
	const float c7 = -2.53113553113553100000000000000000000000E-01f;

	float r = 0.0f, y, w = x;

	do {
		r = r + (1.0f / mc_raise2f(w));
		w = w + 1.0f;
	} while (w < 10.0f);

	y = mc_raise2f(w);
	r = r + (5.0f / w);

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

MC_TARGET_PROC double mc_trigamma_approx1(const double x)
{
//!# Returns psi1(x) by Taylor series expansion.
	const double c1 = +1.0000000000000000000000000000000000000000E+00;
	const double c2 = +1.6666666666666666000000000000000000000000E-01;
	const double c3 = -3.3333333333333333000000000000000000000000E-02;
	const double c4 = +2.3809523809523808000000000000000000000000E-02;
	const double c5 = -3.3333333333333333000000000000000000000000E-02;
	const double c6 = +7.5757575757575760000000000000000000000000E-02;
	const double c7 = -2.5311355311355310000000000000000000000000E-01;

	double r = 0.0, y, w = x;

	do {
		r = r + (1.0 / mc_raise2(w));
		w = w + 1.0;
	} while (w < 10.0);

	y = mc_raise2(w);
	r = r + (5.0 / w);

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

MC_TARGET_PROC long double mc_trigammal_approx1(const long double x)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
//!# Returns psi1(x) by Taylor series expansion.
	const long double c1 = +1.000000000000000000000000000000000000000000000000000000000000000E+00L;
	const long double c2 = +1.666666666666666600000000000000000000000000000000000000000000000E-01L;
	const long double c3 = -3.333333333333333300000000000000000000000000000000000000000000000E-02L;
	const long double c4 = +2.380952380952380800000000000000000000000000000000000000000000000E-02L;
	const long double c5 = -3.333333333333333300000000000000000000000000000000000000000000000E-02L;
	const long double c6 = +7.575757575757576000000000000000000000000000000000000000000000000E-02L;
	const long double c7 = -2.531135531135531000000000000000000000000000000000000000000000000E-01L;

	long double r = 0.0L, y, w = x;

	do {
		r = r + (1.0L / mc_raise2l(w));
		w = w + 1.0L;
	} while (w < 10.0L);

	y = mc_raise2l(w);
	r = r + (5.0L / w);

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
	return mc_cast(long double, mc_trigamma_approx1(mc_cast(const double, x)));
#	endif
}

#pragma mark - mc_trigamma -

MC_TARGET_FUNC float mc_trigammaf(const float x)
{
	return mc_trigammaf_approx1(x);
}

MC_TARGET_FUNC double mc_trigamma(const double x)
{
	return mc_trigamma_approx1(x);
}

MC_TARGET_FUNC long double mc_trigammal(const long double x)
{
	return mc_trigammal_approx1(x);
}

#endif /* !MC_TRIGAMMA_H */

/* EOF */