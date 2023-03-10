//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rgamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_fisint.h>
#include <macadam/details/math/mc_fisneg.h>
#include <macadam/details/math/mc_gammaln.h>
#include <macadam/details/math/mc_gammasign.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_xpolyevaln.h>

#ifndef MC_RGAMMA_H
#define MC_RGAMMA_H

#pragma mark - mc_rgamma_approx0 -

MC_TARGET_FUNC float mc_rgammaf(const float x)
{
//!# Computes reciprocical gamma(x) reducing error rate for low input i.e 1 / gamma(x).
	if (mc_isnan(x) || mc_isinf(x)) {
		if (mc_isinfp(x)) {
			return 0.0f;
		}
		return MCK_NAN;
	}
	if (x == 0.0f) {
		if (mc_fisnegf(x)) {
			return MCK_INFN;
		}
		return MCK_INFP;
	}
	if (x < 0.0f && mc_fisintf(x)) {
		return 0.0f;
	}
	if (x > 0.0f && x < 0.03125f) {
		const float p1 = +1.00000000000000000000000000000000000000E+00f;
		const float p2 = +5.77215664901532860825300000000000000000E-01f;
		const float p3 = -6.55878071520254068466800000000000000000E-01f;
		const float p4 = -4.20026350340334405447300000000000000000E-02f;
		const float p5 = +1.66538611372080520675800000000000000000E-01f;
		const float p6 = -4.21977336070591547008900000000000000000E-02f;
		const float p7 = -9.62202336040627164574400000000000000000E-03f;
		const float p8 = +7.22059947803690967233100000000000000000E-03f;
		const float p9 = -1.19394505138151009561400000000000000000E-03f;

		return x * mc_xpolyeval9f(x, p1, p2, p3, p4, p5, p6, p7, p8, p9);
	}
	if (x < 0.0f && x > -0.03125f) {
		const float p1 = -1.00000000000000000000000000000000000000E+00f;
		const float p2 = +5.77215664901532860872700000000000000000E-01f;
		const float p3 = +6.55878071520253654711600000000000000000E-01f;
		const float p4 = -4.20026350340211291050400000000000000000E-02f;
		const float p5 = -1.66538611394441351933500000000000000000E-01f;
		const float p6 = -4.21977334373119172166400000000000000000E-02f;
		const float p7 = +9.62191115503597673370600000000000000000E-03f;
		const float p8 = +7.22083726189317032570400000000000000000E-03f;
		const float p9 = +1.13337416724389438201000000000000000000E-03f;

		return mc_xpolyeval9f(-x, p1, p2, p3, p4, p5, p6, p7, p8, p9);
	}
	const float g = mc_gammalnf_approx0(x, MC_NULLPTR);
	if (mc_isnan(g) || mc_isinf(g)) {
		return g;
	}
	return mc_gammasignf(x) * mc_expf(-g);
}

MC_TARGET_FUNC double mc_rgamma(const double x)
{
//!# Computes reciprocical gamma(x) reducing error rate for low input i.e 1 / gamma(x).
	if (mc_isnan(x) || mc_isinf(x)) {
		if (mc_isinfp(x)) {
			return 0.0;
		}
		return MCK_NAN;
	}
	if (x == 0.0) {
		if (mc_fisneg(x)) {
			return MCK_INFN;
		}
		return MCK_INFP;
	}
	if (x < 0.0 && mc_fisint(x)) {
		return 0.0;
	}
	if (x > 0.0 && x < 0.03125) {
		const double p1 = +1.0000000000000000000000000000000000000000E+00;
		const double p2 = +5.7721566490153286082530000000000000000000E-01;
		const double p3 = -6.5587807152025406846680000000000000000000E-01;
		const double p4 = -4.2002635034033440544730000000000000000000E-02;
		const double p5 = +1.6653861137208052067580000000000000000000E-01;
		const double p6 = -4.2197733607059154700890000000000000000000E-02;
		const double p7 = -9.6220233604062716457440000000000000000000E-03;
		const double p8 = +7.2205994780369096723310000000000000000000E-03;
		const double p9 = -1.1939450513815100956140000000000000000000E-03;

		return x * mc_xpolyeval9(x, p1, p2, p3, p4, p5, p6, p7, p8, p9);
	}
	if (x < 0.0 && x > -0.03125) {
		const double p1 = -1.0000000000000000000000000000000000000000E+00;
		const double p2 = +5.7721566490153286087270000000000000000000E-01;
		const double p3 = +6.5587807152025365471160000000000000000000E-01;
		const double p4 = -4.2002635034021129105040000000000000000000E-02;
		const double p5 = -1.6653861139444135193350000000000000000000E-01;
		const double p6 = -4.2197733437311917216640000000000000000000E-02;
		const double p7 = +9.6219111550359767337060000000000000000000E-03;
		const double p8 = +7.2208372618931703257040000000000000000000E-03;
		const double p9 = +1.1333741672438943820100000000000000000000E-03;

		return mc_xpolyeval9(-x, p1, p2, p3, p4, p5, p6, p7, p8, p9);
	}
	const double g = mc_gammaln_approx0(x, MC_NULLPTR);
	if (mc_isnan(g) || mc_isinf(g)) {
		return g;
	}
	return mc_gammasign(x) * mc_exp(-g);
}

MC_TARGET_FUNC long double mc_rgammal(const long double x)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
//!# Computes reciprocical gamma(x) reducing error rate for low input i.e 1 / gamma(x).
	if (mc_isnan(x) || mc_isinf(x)) {
		if (mc_isinfp(x)) {
			return 0.0L;
		}
		return MCK_NAN;
	}
	if (x == 0.0L) {
		if (mc_fisnegl(x)) {
			return MCK_INFN;
		}
		return MCK_INFP;
	}
	if (x < 0.0L && mc_fisintl(x)) {
		return 0.0L;
	}
	if (x > 0.0L && x < 0.03125L) {
		const long double p1 = +1.000000000000000000000000000000000000000000000000000000000000000E+00L;
		const long double p2 = +5.772156649015328608253000000000000000000000000000000000000000000E-01L;
		const long double p3 = -6.558780715202540684668000000000000000000000000000000000000000000E-01L;
		const long double p4 = -4.200263503403344054473000000000000000000000000000000000000000000E-02L;
		const long double p5 = +1.665386113720805206758000000000000000000000000000000000000000000E-01L;
		const long double p6 = -4.219773360705915470089000000000000000000000000000000000000000000E-02L;
		const long double p7 = -9.622023360406271645744000000000000000000000000000000000000000000E-03L;
		const long double p8 = +7.220599478036909672331000000000000000000000000000000000000000000E-03L;
		const long double p9 = -1.193945051381510095614000000000000000000000000000000000000000000E-03L;

		return x * mc_xpolyeval9l(x, p1, p2, p3, p4, p5, p6, p7, p8, p9);
	}
	if (x < 0.0L && x > -0.03125L) {
		const long double p1 = -1.000000000000000000000000000000000000000000000000000000000000000E+00L;
		const long double p2 = +5.772156649015328608727000000000000000000000000000000000000000000E-01L;
		const long double p3 = +6.558780715202536547116000000000000000000000000000000000000000000E-01L;
		const long double p4 = -4.200263503402112910504000000000000000000000000000000000000000000E-02L;
		const long double p5 = -1.665386113944413519335000000000000000000000000000000000000000000E-01L;
		const long double p6 = -4.219773343731191721664000000000000000000000000000000000000000000E-02L;
		const long double p7 = +9.621911155035976733706000000000000000000000000000000000000000000E-03L;
		const long double p8 = +7.220837261893170325704000000000000000000000000000000000000000000E-03L;
		const long double p9 = +1.133374167243894382010000000000000000000000000000000000000000000E-03L;

		return mc_xpolyeval9l(-x, p1, p2, p3, p4, p5, p6, p7, p8, p9);
	}
	const long double g = mc_gammalnl_approx0(x, MC_NULLPTR);
	if (mc_isnan(g) || mc_isinf(g)) {
		return g;
	}
	return mc_gammasignl(x) * mc_expl(-g);
#	else
	return mc_cast(long double, mc_rgamma(mc_cast(double, x)));
#	endif
}

#endif /* !MC_RGAMMA_H */

/* EOF */