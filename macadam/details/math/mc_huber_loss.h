//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_huber_loss.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_rsqrt.h>

#ifndef MC_HUBER_LOSS_H
#define MC_HUBER_LOSS_H

#pragma mark - mc_huber_loss -

MC_TARGET_FUNC float mc_huber_lossf(const float r, const float c, const float s, const int p)
{
//!# Huber loss functions. Pass p=1 for pseudo.
//!# Default settings c=1.345 and s=1 (scale)
	float d;
	if (p) {
		const float w = s != 0.0f ? 1.0f / s : 1.0f;
		const float h = 1.0f / (w * c);
		d             = mc_rsqrtf(1.0f + mc_raise2f(r * h));
	} else {
		const float w = s != 0.0f ? 1.0f / s : 1.0f;
		const float h = mc_fabsf(r * w);
		d             = h < c ? 1.0f : c * (1.0f / h);
	}
	return d;
}

MC_TARGET_FUNC double mc_huber_loss(const double r, const double c, const double s, const int p)
{
//!# Huber loss functions. Pass p=1 for pseudo.
//!# Default settings c=1.345 and s=1 (scale)
	double d;
	if (p) {
		const double w = s != 0.0 ? 1.0 / s : 1.0;
		const double h = 1.0 / (w * c);
		d              = mc_rsqrt(1.0 + mc_raise2(r * h));
	} else {
		const double w = s != 0.0 ? 1.0 / s : 1.0;
		const double h = mc_fabs(r * w);
		d              = h < c ? 1.0 : c * (1.0 / h);
	}
	return d;
}

MC_TARGET_FUNC long double mc_huber_lossl(const long double r, const long double c, const long double s, const int p)
{
//!# Huber loss functions. Pass p=1 for pseudo.
//!# Default settings c=1.345 and s=1 (scale)
	long double d;
	if (p) {
		const long double w = s != 0.0L ? 1.0L / s : 1.0L;
		const long double h = 1.0L / (w * c);
		d                   = mc_rsqrtl(1.0L + mc_raise2l(r * h));
	} else {
		const long double w = s != 0.0L ? 1.0L / s : 1.0L;
		const long double h = mc_fabsl(r * w);
		d                   = h < c ? 1.0L : c * (1.0L / h);
	}
	return d;
}

#endif /* !MC_HUBER_LOSS_H */

/* EOF */