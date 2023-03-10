//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ibeta.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_lbeta.h>
#include <macadam/details/math/mc_log.h>

#ifndef MC_RBETA_H
#define MC_RBETA_H

#pragma mark - mc_ibeta -

/*! Computes the regularized incomplete beta function.
 *
 * \brief Computing the regularized incomplete beta function.
 *
 * float mc_ibetaf(const float a, const float b, const float x);
 *
 * \param a Value strictly positive a>0.
 * \param b Value strictly positive b>0.
 * \param x Value in the range [0,1].
 * \result  The normalised incomplete beta function of a, b and x.
 */
MC_TARGET_FUNC float mc_ibetaf(const float a, const float b, const float x)
{
	unsigned int i = 0, j;
	float g, k, w, f = 1.0f, c = 1.0f, d = 0.0f;
//!# NAN input results in NAN output.
	if (mc_isnan(a) || mc_isnan(b) || mc_isnan(x)) {
		return MCK_NAN;
	}
//!# x is out of bounds hence infinity is returned.
	if (x < 0.0f || x > 1.0f) {
		return x < 0.0f ? MCK_INFN : MCK_INFP;
	}
	if (x > (a + 1.0f) / (a + b + 2.0f)) {
	//!# The beta inverse is symetric.
		return (1.0f - mc_ibetaf(b, a, 1.0f - x));
	}
//!# Computing delta-gamma + front integral.
	g = mc_lbetaf(a, b);
	k = mc_expf(mc_logf(x) * a + mc_logf(1.0f - x) * b - g) / a;
//!# Reducing, converging.
	for (; i < 256; i++) {
		j = i / 2;
		if (i == 0) {
		//!# First iteration.
			w = 1.0f;
		} else if ((i % 2) == 0) {
			w = (j * (b - j) * x) / ((a + 2.0f * j - 1.0f) * (a + 2.0f * j));
		} else {
			w = -((a + j) * (a + b + j) * x) / ((a + 2.0f * j) * (a + 2.0f * j + 1));
		}
		d = 1.0f + w * d;
		if (mc_fabsf(d) < 1.0E-7f) {
		//!# Clipping to absolute min.
			d = 1.0E-7f;
		}
		d = 1.0f / d;
		c = 1.0f + w / c;
		if (mc_fabsf(c) < 1.0E-7f) {
		//!# Clipping to absolute min.
			c = 1.0E-7f;
		}
		f = f * (c * d);
		if (mc_fabsf(1.0f - (c * d)) < 1.0E-5f) {
			return k * (f - 1.0f);
		}
	}
//!# Unable to reduce, returning towards infinity.
	return MCK_INFP;
}

/*! Computes the regularized incomplete beta function.
 *
 * \brief Computing the regularized incomplete beta function.
 *
 * double mc_ibeta(const double a, const double b, const double x);
 *
 * \param a Value strictly positive a>0.
 * \param b Value strictly positive b>0.
 * \param x Value in the range [0,1].
 * \result  The normalised incomplete beta function of a, b and x.
 */
MC_TARGET_FUNC double mc_ibeta(const double a, const double b, const double x)
{
	unsigned int i = 0, j;
	double g, k, w, f = 1.0, c = 1.0, d = 0.0;
//!# NAN input results in NAN output.
	if (mc_isnan(a) || mc_isnan(b) || mc_isnan(x)) {
		return MCK_NAN;
	}
//!# x is out of bounds hence infinity is returned.
	if (x < 0.0 || x > 1.0) {
		return x < 0.0 ? MCK_INFN : MCK_INFP;
	}
	if (x > (a + 1.0) / (a + b + 2.0)) {
	//!# The beta inverse is symetric.
		return (1.0 - mc_ibeta(b, a, 1.0 - x));
	}
//!# Computing delta-gamma + front integral.
	g = mc_lbeta(a, b);
	k = mc_exp(mc_log(x) * a + mc_log(1.0 - x) * b - g) / a;
//!# Reducing, converging.
	for (; i < 256; i++) {
		j = i / 2;
		if (i == 0) {
		//!# First iteration.
			w = 1.0;
		} else if ((i % 2) == 0) {
			w = (j * (b - j) * x) / ((a + 2.0 * j - 1.0) * (a + 2.0 * j));
		} else {
			w = -((a + j) * (a + b + j) * x) / ((a + 2.0 * j) * (a + 2.0 * j + 1));
		}
		d = 1.0 + w * d;
		if (mc_fabs(d) < 1.0E-15) {
		//!# Clipping to absolute min.
			d = 1.0E-15;
		}
		d = 1.0 / d;
		c = 1.0 + w / c;
		if (mc_fabs(c) < 1.0E-15) {
		//!# Clipping to absolute min.
			c = 1.0E-15;
		}
		f = f * (c * d);
		if (mc_fabs(1.0 - (c * d)) < 1.0E-8) {
			return k * (f - 1.0);
		}
	}
//!# Unable to reduce, returning towards infinity.
	return MCK_INFP;
}

/*! Computes the regularized incomplete beta function.
 *
 * \brief Computing the regularized incomplete beta function.
 *
 * long double mc_ibetal(const long double a, const long double b, const long double x);
 *
 * \param a Value strictly positive a>0.
 * \param b Value strictly positive b>0.
 * \param x Value in the range [0,1].
 * \result  The normalised incomplete beta function of a, b and x.
 */
MC_TARGET_FUNC long double mc_ibetal(const long double a, const long double b, const long double x)
{
	unsigned int i = 0, j;
	long double g, k, w, f = 1.0L, c = 1.0L, d = 0.0L;
//!# NAN input results in NAN output.
	if (mc_isnan(a) || mc_isnan(b) || mc_isnan(x)) {
		return MCK_NAN;
	}
//!# x is out of bounds hence infinity is returned.
	if (x < 0.0L || x > 1.0L) {
		return x < 0.0L ? MCK_INFN : MCK_INFP;
	}
	if (x > (a + 1.0L) / (a + b + 2.0L)) {
	//!# The beta inverse is symetric.
		return (1.0L - mc_ibetal(b, a, 1.0L - x));
	}
//!# Computing delta-gamma + front integral.
	g = mc_lbetal(a, b);
	k = mc_expl(mc_logl(x) * a + mc_logl(1.0L - x) * b - g) / a;
//!# Reducing, converging.
	for (; i < 256; i++) {
		j = i / 2;
		if (i == 0) {
		//!# First iteration.
			w = 1.0L;
		} else if ((i % 2) == 0) {
			w = (j * (b - j) * x) / ((a + 2.0L * j - 1.0L) * (a + 2.0L * j));
		} else {
			w = -((a + j) * (a + b + j) * x) / ((a + 2.0L * j) * (a + 2.0L * j + 1));
		}
		d = 1.0L + w * d;
		if (mc_fabsl(d) < 1.0E-30L) {
		//!# Clipping to absolute min.
			d = 1.0E-30L;
		}
		d = 1.0L / d;
		c = 1.0L + w / c;
		if (mc_fabsl(c) < 1.0E-30L) {
		//!# Clipping to absolute min.
			c = 1.0E-30L;
		}
		f = f * (c * d);
		if (mc_fabsl(1.0L - (c * d)) < 1.0E-12L) {
			return k * (f - 1.0L);
		}
	}
//!# Unable to reduce, returning towards infinity.
	return MCK_INFP;
}

#endif /* !MC_RBETA_H */

/* EOF */