//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_inverf.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_erf.h>
#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_INVERF_H
#define MC_INVERF_H

#pragma mark - mc_inverf -

MC_TARGET_FUNC float mc_inverff(const float x)
{
	float w = 0.0f, z, n, d;
	if (x < -1.0f || x > 1.0f) {
		return MCK_NAN;
	} else if (x == -1.0f) {
		return MCK_INFP;
	} else if (x == 1.0f) {
		return MCK_INFN;
	}
	if (mc_fabsf(x) <= 0.7f) {
		z = mc_raise2f(x);
		n = ((((-0.140543331f * z + 0.914624893f) * z - 1.645349621f) * z + 0.886226899f));
		d = (((( 0.012229801f * z - 0.329097515f) * z + 1.442710462f) * z - 2.118377725f) * z + 1.0f);
		w = (x * n) / d;
	} else if ((mc_fabsf(x) > 0.7f) && (mc_fabsf(x) < 1.0f) ) {
		z = mc_sqrtf(-mc_logf((1.0f - mc_fabsf(x)) * 0.5f));
		n = ((1.641345311f * z + 3.429567803f) * z - 1.624906493f) * z - 1.970840454f;
		d = ((1.637067800f * z + 3.543889200f) * z + 1.0f);
		w = (x < 0 ? -n : n) / d;
	}
	w = w - ((mc_erff(w) - x) / (MCK_KF(MCK_2_SQRTPI) * mc_expf(-w * w)));
	w = w - ((mc_erff(w) - x) / (MCK_KF(MCK_2_SQRTPI) * mc_expf(-w * w)));
	return w;
}

MC_TARGET_FUNC double mc_inverf(const double x)
{
	double w = 0.0, z, n, d;
	if (x < -1.0 || x > 1.0) {
		return MCK_NAN;
	} else if (x == -1.0) {
		return MCK_INFP;
	} else if (x == 1.0) {
		return MCK_INFN;
	}
	if (mc_fabs(x) <= 0.7) {
		z = mc_raise2(x);
		n = ((((-0.140543331 * z + 0.914624893) * z - 1.645349621) * z + 0.886226899));
		d = (((( 0.012229801 * z - 0.329097515) * z + 1.442710462) * z - 2.118377725) * z + 1.0);
		w = (x * n) / d;
	} else if ((mc_fabs(x) > 0.7) && (mc_fabs(x) < 1.0) ) {
		z = mc_sqrt(-mc_log((1.0 - mc_fabs(x)) * 0.5));
		n = ((1.641345311 * z + 3.429567803) * z - 1.624906493) * z - 1.970840454;
		d = ((1.637067800 * z + 3.543889200) * z + 1.0);
		w = (x < 0 ? -n : n) / d;
	}
	w = w - ((mc_erf(w) - x) / (MCK_K(MCK_2_SQRTPI) * mc_exp(-w * w)));
	w = w - ((mc_erf(w) - x) / (MCK_K(MCK_2_SQRTPI) * mc_exp(-w * w)));
	return w;
}

MC_TARGET_FUNC long double mc_inverfl(const long double x)
{
	long double w = 0.0L, z, n, d;
	if (x < -1.0L || x > 1.0L) {
		return MCK_NAN;
	} else if (x == -1.0L) {
		return MCK_INFP;
	} else if (x == 1.0L) {
		return MCK_INFN;
	}
	if (mc_fabsl(x) <= 0.7L) {
		z = mc_raise2l(x);
		n = ((((-0.140543331L * z + 0.914624893L) * z - 1.645349621L) * z + 0.886226899L));
		d = (((( 0.012229801L * z - 0.329097515L) * z + 1.442710462L) * z - 2.118377725L) * z + 1.0L);
		w = (x * n) / d;
	} else if ((mc_fabsl(x) > 0.7L) && (mc_fabsl(x) < 1.0L) ) {
		z = mc_sqrtl(-mc_logl((1.0L - mc_fabsl(x)) * 0.5L));
		n = ((1.641345311L * z + 3.429567803L) * z - 1.624906493L) * z - 1.970840454L;
		d = ((1.637067800L * z + 3.543889200L) * z + 1.0L);
		w = (x < 0 ? -n : n) / d;
	}
	w = w - ((mc_erfl(w) - x) / (MCK_KL(MCK_2_SQRTPI) * mc_expl(-w * w)));
	w = w - ((mc_erfl(w) - x) / (MCK_KL(MCK_2_SQRTPI) * mc_expl(-w * w)));
	return w;
}

#endif /* !MC_INVERF_H */

/* EOF */