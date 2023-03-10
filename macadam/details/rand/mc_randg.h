//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_randg.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zsqrt.h>
#include <macadam/details/rand/mc_rand_uniform.h>

#ifndef MC_RANDG_H
#define MC_RANDG_H

#pragma mark - mc_randg -

MC_TARGET_PROC float mc_randgf(void)
{
#	if MC_TARGET_RAND_USE_BOXMULLER
//!# Box-Muller transform. Standard gaussian (normal) distribution for mean=0, stddev=1.
	static volatile int phase_s = 0;
	static volatile float x_s   = 0.0f;
	float r = 0.0f, u, v;

	if (phase_s != 0) {
		r = x_s;
	} else {
		do {
			const float r1 = mc_randuf();
			const float r2 = mc_randuf();
			u              = r1;
			v              = r2;
		}
		while (u <= MCLIMITS_EPSILONF * 2.0f);
		r   = mc_sqrtf(-2.0f * mc_logf(u)) * mc_cosf(MCK_KF(MCK_2PI) * v);
		x_s = mc_sqrtf(-2.0f * mc_logf(u)) * mc_sinf(MCK_KF(MCK_2PI) * v);
	}
	phase_s = !phase_s;
	return r;
#	else
//!# Marsaglia polar transform. Standard gaussian (normal) distribution for mean=0, stddev=1.
	static volatile int phase_s = 0;
	static volatile float x_s   = 0.0f;
	float r, s = 0.0f, u, v, w;
	if (phase_s != 0) {
		r = x_s;
	} else {
		do {
			const float r1 = mc_randuf();
			const float r2 = mc_randuf();
			u              = 2.0f * r1 - 1.0f;
			v              = 2.0f * r2 - 1.0f;
			s              = mc_raise2f(u) + mc_raise2f(v);
		} while (s > 1.0f || s == 0.0f);
		w = -2.0f * mc_logf(s) * (1.0f / s);
		if (w < 0.0f) {
			mc_zsqrtf(&w, &s, w, 0.0f);
		} else {
			s = mc_sqrtf(w);
		}
		r   = u * s;
		x_s = v * s;
	}
	phase_s = !phase_s;
	return r;
#	endif
}

MC_TARGET_PROC double mc_randg(void)
{
#	if MC_TARGET_RAND_USE_BOXMULLER
//!# Box-Muller transform. Standard gaussian (normal) distribution for mean=0, stddev=1.
	static volatile int phase_s = 0;
	static volatile double x_s  = 0.0;
	double r = 0.0, u, v;

	if (phase_s != 0) {
		r = x_s;
	} else {
		do {
			const double r1 = mc_randu();
			const double r2 = mc_randu();
			u               = r1;
			v               = r2;
		}
		while (u <= MCLIMITS_EPSILON * 2.0);
		r   = mc_sqrt(-2.0 * mc_log(u)) * mc_cos(MCK_K(MCK_2PI) * v);
		x_s = mc_sqrt(-2.0 * mc_log(u)) * mc_sin(MCK_K(MCK_2PI) * v);
	}
	phase_s = !phase_s;
	return r;
#	else
//!# Marsaglia polar transform. Standard gaussian (normal) distribution for mean=0, stddev=1.
	static volatile int phase_s = 0;
	static volatile double x_s  = 0.0;
	double r, s = 0.0, u, v, w;
	if (phase_s != 0) {
		r = x_s;
	} else {
		do {
			const double r1 = mc_randu();
			const double r2 = mc_randu();
			u               = 2.0 * r1 - 1.0;
			v               = 2.0 * r2 - 1.0;
			s               = mc_raise2(u) + mc_raise2(v);
		} while (s > 1.0 || s == 0.0);
		w = -2.0 * mc_log(s) * (1.0 / s);
		if (w < 0.0) {
			mc_zsqrt(&w, &s, w, 0.0);
		} else {
			s = mc_sqrt(w);
		}
		r   = u * s;
		x_s = v * s;
	}
	phase_s = !phase_s;
	return r;
#	endif
}

MC_TARGET_PROC long double mc_randgl(void)
{
#	if MC_TARGET_RAND_USE_BOXMULLER
//!# Box-Muller transform. Standard gaussian (normal) distribution for mean=0, stddev=1.
	static volatile int phase_s     = 0;
	static volatile long double x_s = 0.0L;
	long double r = 0.0L, u, v;
s
	if (phase_s != 0) {
		r = x_s;
	} else {
		do {
			const long double r1 = mc_randul();
			const long double r2 = mc_randul();
			u                    = r1;
			v                    = r2;
		}
		while (u <= MCLIMITS_EPSILONL  * 2.0L);
		r   = mc_sqrtl(-2.0L * mc_logl(u)) * mc_cosl(MCK_KL(MCK_2PI) * v);
		x_s = mc_sqrtl(-2.0L * mc_logl(u)) * mc_sinl(MCK_KL(MCK_2PI) * v);
	}
	phase_s = !phase_s;
	return r;
#	else
//!# Marsaglia polar transform. Standard gaussian (normal) distribution for mean=0, stddev=1.
	static volatile int phase_s     = 0;
	static volatile long double x_s = 0.0L;
	long double r, s = 0.0L, u, v, w;
	if (phase_s != 0) {
		r = x_s;
	} else {
		do {
			const long double r1 = mc_randul();
			const long double r2 = mc_randul();
			u                    = 2.0L * r1 - 1.0L;
			v                    = 2.0L * r2 - 1.0L;
			s                    = mc_raise2l(u) + mc_raise2l(v);
		} while (s > 1.0L || s == 0.0L);
		w = -2.0L * mc_logl(s) * (1.0L / s);
		if (w < 0.0L) {
			mc_zsqrtl(&w, &s, w, 0.0L);
		} else {
			s = mc_sqrtl(w);
		}
		r   = u * s;
		x_s = v * s;
	}
	phase_s = !phase_s;
	return r;
#	endif
}

#endif /* !MC_RANDG_H */

/* EOF */