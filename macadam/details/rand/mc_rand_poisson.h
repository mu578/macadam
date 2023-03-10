//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_poisson.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_floor.h>
#include <macadam/details/math/mc_gammaln.h>
#include <macadam/details/math/mc_itrunc32.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_tanpi.h>
#include <macadam/details/rand/mc_randu.h>

#ifndef MC_RAND_POISSON_H
#define MC_RAND_POISSON_H

#pragma mark - mc_rand_poisson -

MC_TARGET_FUNC int rand_poissonf(const float l)
{
	float r, m = 0.0f, p = 1.0f, x, f;

	if (l < 90.0f) {
		r = mc_expf(-l);
		do {
			m = m + 1.0f;
			p = p * mc_randuf();
		} while (p > r);
		m = m - 1.0f;
	} else {
//!# rejection can be very slow.
		const float w = mc_sqrtf(l);
		const float g = mc_logf(l);
		x             = -1.0f;
		do {
			do {
				p = mc_randuf();
				x = l + w * mc_tanpif(p - 1.0f * 0.5f);
			} while (x < 0.0f);
			m = mc_floorf(x);
			p = w / (MCK_KF(MCK_PI) * ((x - l) * (x - l) + l));
			f = m * g - l - mc_gammalnf(m + 1.0f);
			f = mc_expf(f);
			r = f / p / 2.4f;
		} while (mc_randuf() > r);
	}
	return mc_itrunc32f(m);
}

MC_TARGET_FUNC int rand_poisson(const double l)
{
	double r, m = 0.0, p = 1.0, x, f;

	if (l < 180.0) {
		r = mc_exp(-l);
		do {
			m = m + 1.0;
			p = p * mc_randu();
		} while (p > r);
		m = m - 1.0;
	} else {
//!# rejection can be very slow.
		const double w = mc_sqrt(l);
		const double g = mc_log(l);
		x              = -1.0;
		do {
			do {
				p = mc_randu();
				x = l + w * mc_tanpi(p - 1.0 * 0.5);
			} while (x < 0.0);
			m = mc_floor(x);
			p = w / (MCK_K(MCK_PI) * ((x - l) * (x - l) + l));
			f = m * g - l - mc_gammaln(m + 1.0);
			f = mc_exp(f);
			r = f / p / 2.4;
		} while (mc_randu() > r);
	}
	return mc_itrunc32(m);
}

MC_TARGET_FUNC int rand_poissonl(const long double l)
{
	long double r, m = 0.0L, p = 1.0L, x, f;

	if (l < 720.0L) {
		r = mc_expl(-l);
		do {
			m = m + 1.0L;
			p = p * mc_randul();
		} while (p > r);
		m = m - 1.0;
	} else {
//!# rejection can be very slow.
		const long double w = mc_sqrtl(l);
		const long double g = mc_logl(l);
		x                   = -1.0L;
		do {
			do {
				p = mc_randul();
				x = l + w * mc_tanpil(p - 1.0L * 0.5L);
			} while (x < 0.0L);
			m = mc_floorl(x);
			p = w / (MCK_KL(MCK_PI) * ((x - l) * (x - l) + l));
			f = m * g - l - mc_gammalnl(m + 1.0L);
			f = mc_expl(f);
			r = f / p / 2.4L;
		} while (mc_randul() > r);
	}
	return mc_itrunc32l(m);
}

#endif /* !MC_RAND_POISSON_H */

/* EOF */