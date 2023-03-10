//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_circular.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_acos.h>
#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/rand/mc_randg.h>
#include <macadam/details/rand/mc_randu.h>

#ifndef MC_RAND_CIRCULAR_H
#define MC_RAND_CIRCULAR_H

#pragma mark - mc_rand_circular -

MC_TARGET_FUNC float mc_rand_circularf(const float mu, const float k)
{
//!# Circular normal/uniform distribution generator.
	float x = MCK_KF(MCK_PI);
	if (k >= 0.0f) {
		if (k == 0.0f) {
//!# Circular uniform.
			x = mc_randuf() * MCK_KF(MCK_2PI);
		} else if (mc_fabsf(mu) <= MCK_KF(MCK_PI)) {
//!# Circular normal.
			const float a = (1.0f + mc_sqrtf(1.0f + 4.0f * mc_raise2f(k))) / (2.0f * k);
			for (;;) {
				const float u = mc_randgf();
				const float v = mc_randgf();
				const float g = MCK_KF(MCK_2PI) * u;
				const float z = mc_cosf(g);
				const float w = (1.0f - a * z) / (a - z);
				const float t = k * (a - w);
				const float e = mc_expf(1.0f - t);
				if (!(v > t * (2.0f - t) && v > t * e)) {
					const float y = u < 0.5f ? -mc_acosf(w) : mc_acosf(w);
					x             = y + mu;
					x             = x >= MCK_KF(MCK_PI) ? x - MCK_KF(MCK_PI) : (x < -MCK_KF(MCK_PI) ? x + MCK_KF(MCK_PI) : x);
					break;
				}
			}
		}
	}
	return x;
}

MC_TARGET_FUNC double mc_rand_circularff(const float mu, const float k)
{
//!# Circular normal/uniform distribution generator.
	const double mud = mc_cast(const double, mu);
	const double kd  = mc_cast(const double, k);
	double x         = MCK_K(MCK_PI);
	if (kd >= 0.0) {
		if (kd == 0.0) {
//!# Circular uniform.
			x = mc_randu() * MCK_K(MCK_2PI);
		} else if (mc_fabs(mud) <= MCK_K(MCK_PI)) {
//!# Circular normal.
			const double a = (1.0 + mc_sqrt(1.0 + 4.0 * mc_raise2(kd))) / (2.0 * kd);
			for (;;) {
				const double u = mc_randg();
				const double v = mc_randg();
				const double g = MCK_K(MCK_2PI) * u;
				const double z = mc_cos(g);
				const double w = (1.0 - a * z) / (a - z);
				const double t = k * (a - w);
				const double e = mc_exp(1.0 - t);
				if (!(v > t * (2.0 - t) && v > t * e)) {
					const double y = u < 0.5 ? -mc_acos(w) : mc_acos(w);
					x              = y + mud;
					x              = x >= MCK_K(MCK_PI) ? x - MCK_K(MCK_PI) : (x < -MCK_K(MCK_PI) ? x + MCK_K(MCK_PI) : x);
					break;
				}
			}
		}
	}
	return x;
}

MC_TARGET_FUNC double mc_rand_circular(const double mu, const double k)
{
//!# Circular normal/uniform distribution generator.
	double x = MCK_K(MCK_PI);
	if (k >= 0.0) {
		if (k == 0.0) {
//!# Circular uniform.
			x = mc_randu() * MCK_K(MCK_2PI);
		} else if (mc_fabs(mu) <= MCK_K(MCK_PI)) {
//!# Circular normal.
			const double a = (1.0 + mc_sqrt(1.0 + 4.0 * mc_raise2(k))) / (2.0 * k);
			for (;;) {
				const double u = mc_randg();
				const double v = mc_randg();
				const double g = MCK_K(MCK_2PI) * u;
				const double z = mc_cos(g);
				const double w = (1.0 - a * z) / (a - z);
				const double t = k * (a - w);
				const double e = mc_exp(1.0 - t);
				if (!(v > t * (2.0 - t) && v > t * e)) {
					const double y = u < 0.5 ? -mc_acos(w) : mc_acos(w);
					x              = y + mu;
					x              = x >= MCK_K(MCK_PI) ? x - MCK_K(MCK_PI) : (x < -MCK_K(MCK_PI) ? x + MCK_K(MCK_PI) : x);
					break;
				}
			}
		}
	}
	return x;
}

MC_TARGET_FUNC long double mc_rand_circularl(const long double mu, const long double k)
{
//!# Circular normal/uniform distribution generator.
	long double x = MCK_KL(MCK_PI);
	if (k >= 0.0L) {
		if (k == 0.0L) {
//!# Circular uniform.
			x = mc_randul() * MCK_KL(MCK_2PI);
		} else if (mc_fabsl(mu) <= MCK_KL(MCK_PI)) {
//!# Circular normal.
			const long double a = (1.0L + mc_sqrtl(1.0L + 4.0L * mc_raise2l(k))) / (2.0L * k);
			for (;;) {
				const long double u = mc_randgl();
				const long double v = mc_randgl();
				const long double g = MCK_KL(MCK_2PI) * u;
				const long double z = mc_cosl(g);
				const long double w = (1.0L - a * z) / (a - z);
				const long double t = k * (a - w);
				const long double e = mc_expl(1.0L - t);
				if (!(v > t * (2.0L - t) && v > t * e)) {
					const long double y = u < 0.5L ? -mc_acosl(w) : mc_acosl(w);
					x                   = y + mu;
					x                   = x >= MCK_KL(MCK_PI) ? x - MCK_KL(MCK_PI) : (x < -MCK_KL(MCK_PI) ? x + MCK_KL(MCK_PI) : x);
					break;
				}
			}
		}
	}
	return x;
}

#endif /* !MC_RAND_CIRCULAR_H */

/* EOF */