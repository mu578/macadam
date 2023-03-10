//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_lgamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_raise3.h>
#include <macadam/details/math/mc_pow.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/rand/mc_randg.h>
#include <macadam/details/rand/mc_randu.h>

#ifndef MC_RAND_LGAMMA_H
#define MC_RAND_LGAMMA_H

#pragma mark - mc_rand_lgamma -

MC_TARGET_FUNC float mc_rand_lgammaf(const float a, const float l)
{
//!# Log(Gamma(alpha,lambda)) generator using Marsaglia and Tsang method.
//!# Log Gamma RNG a=alpha=shape, l=lambda=scale.
	float r = 0.0f, b, d, c, v, u;
	unsigned int j = 0;
	if (a <= 0.0f || l <= 0.0f) {
		return r;
	}
	b = a;
	if (b < 1.0f) {
		j = 1;
		b = b + 1.0f;
	}
	d = b - 1.0f / 3.0f;
	c = 1.0f / 3.0f / mc_sqrtf(d);
	do {
		do {
			r = mc_randgf();
			v = 1.0f + c * r;
		} while ( v <= 0 );
		v = mc_raise3f(v);
		r = mc_raise2f(r);
		u = mc_randuf();
		if (u <= 1.0f - 0.331f * mc_raise2f(r)) {
			break;
		}
		u = mc_logf(u);
	} while (!(u <= 0.5f * r + d * (1.0f - v + mc_logf(v))));
	r = d * v;
	if (j) {
		r = r * mc_powf(1.0f -  mc_randuf(), 1.0f / (b - 1.0f));
	}
	return mc_logf(r / l);
}

MC_TARGET_FUNC double mc_rand_lgammaff(const float a, const float l)
{
//!# Log(Gamma(alpha,lambda)) generator using Marsaglia and Tsang method.
//!# Log Gamma RNG a=alpha=shape, l=lambda=scale.
	double r = 0.0, b, d, c, v, u;
	unsigned int j = 0;
	if (a <= 0.0f || l <= 0.0f) {
		return r;
	}
	b = mc_cast(double, a);
	if (b < 1.0) {
		j = 1;
		b = b + 1.0;
	}
	d = b - 1.0 / 3.0;
	c = 1.0 / 3.0 / mc_sqrt(d);
	do {
		do {
			r = mc_randg();
			v = 1.0 + c * r;
		} while ( v <= 0 );
		v = mc_raise3(v);
		r = mc_raise2(r);
		u = mc_randu();
		if (u <= 1.0 - 0.331 * mc_raise2(r)) {
			break;
		}
		u = mc_log(u);
	} while (!(u <= 0.5 * r + d * (1.0 - v + mc_log(v))));
	r = d * v;
	if (j) {
		r = r * mc_pow(1.0 -  mc_randu(), 1.0 / (b - 1.0));
	}
	return mc_log(r / mc_cast(const double, l));
}

MC_TARGET_FUNC double mc_rand_lgamma(const double a, const double l)
{
//!# Log(Gamma(alpha,lambda)) generator using Marsaglia and Tsang method.
//!# Log Gamma RNG a=alpha=shape, l=lambda=scale.
	double r = 0.0, b, d, c, v, u;
	unsigned int j = 0;
	if (a <= 0.0 || l <= 0.0) {
		return r;
	}
	b = a;
	if (b < 1.0) {
		j = 1;
		b = b + 1.0;
	}
	d = b - 1.0 / 3.0;
	c = 1.0 / 3.0 / mc_sqrt(d);
	do {
		do {
			r = mc_randg();
			v = 1.0 + c * r;
		} while ( v <= 0 );
		v = mc_raise3(v);
		r = mc_raise2(r);
		u = mc_randu();
		if (u <= 1.0 - 0.331 * mc_raise2(r)) {
			break;
		}
		u = mc_log(u);
	} while (!(u <= 0.5 * r + d * (1.0 - v + mc_log(v))));
	r = d * v;
	if (j) {
		r = r * mc_pow(1.0 -  mc_randu(), 1.0 / (b - 1.0));
	}
	return mc_log(r / l);
}

MC_TARGET_FUNC long double mc_rand_lgammal(const long double a, const long double l)
{
//!# Log(Gamma(alpha,lambda)) generator using Marsaglia and Tsang method.
//!# Log Gamma RNG a=alpha=shape, l=lambda=scale.
	long double r = 0.0L, b, d, c, v, u;
	unsigned int j = 0;
	if (a <= 0.0L || l <= 0.0L) {
		return r;
	}
	b = a;
	if (b < 1.0L) {
		j = 1;
		b = b + 1.0L;
	}
	d = b - 1.0L / 3.0L;
	c = 1.0L / 3.0L / mc_sqrtl(d);
	do {
		do {
			r = mc_randgl();
			v = 1.0L + c * r;
		} while ( v <= 0 );
		v = mc_raise3l(v);
		r = mc_raise2l(r);
		u = mc_randul();
		if (u <= 1.0L - 0.331L * mc_raise2l(r)) {
			break;
		}
		u = mc_logl(u);
	} while (!(u <= 0.5L * r + d * (1.0L - v + mc_logl(v))));
	r = d * v;
	if (j) {
		r = r * mc_powl(1.0L -  mc_randul(), 1.0L / (b - 1.0L));
	}
	return mc_logl(r / l);
}

#endif /* !MC_RAND_LGAMMA_H */

/* EOF */