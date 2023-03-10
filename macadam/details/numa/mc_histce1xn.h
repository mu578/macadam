//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_histce1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_logbase.h>

#ifndef MC_HISTCE1XN_H
#define MC_HISTCE1XN_H

#pragma mark - mc_histce1xn -

MC_TARGET_FUNC float mc_histce1xnf(int npts, int nbins, int logbase, const int * h)
{
//!# Requires h[npts] where 1 < n and 0 < nbins and 0 < b. Computing the entropy
//!# of a given `counted` histogram (frequency) @see  `mc_histcs1xn` or `mc_histcg1xn`.
//!#     npts    - Size of histogram i.e h[npts]. @note npts, nbins are usually equal @see `mc_histcs1xn`.
//!#     nbins   - The given bin width.           @note npts, nbins are usually equal @see `mc_histcs1xn`.
//!#     logbase - The logarithm base, if set to `one`, loge is used i.e nat; `two` for bit unit.

	int i       = 0;
	const int c = nbins < 1 ? 1 : nbins;
	float e     = 0.0f;

	if (npts > 1 && logbase > 0) {
		if (c < 0x1000001) {
			for (; i < npts; i++) {
				const float p  = mc_cast(float, h[i]) / c;
				e             -= p > 0.0f ? p * mc_logbasef(p, logbase) : 0.0f;
			}
		} else {
			for (; i < npts; i++) {
				const double p  = mc_cast(double, h[i]) / c;
				e              -= mc_cast_expr(float, p > 0.0 ? p * mc_logbase(p, logbase) : 0.0);
			}
		}
	}
	return e;
}

MC_TARGET_FUNC double mc_histce1xn(int npts, int nbins, int logbase, const int * h)
{
//!# Requires h[npts] where 1 < n and 0 < nbins and 0 < b. Computing the entropy
//!# of a given `counted` histogram (frequency) @see  `mc_histcs1xn` or `mc_histcg1xn`.
//!#     npts    - Size of histogram i.e h[npts]. @note npts, nbins are usually equal @see `mc_histcs1xn`.
//!#     nbins   - The given bin width.           @note npts, nbins are usually equal @see `mc_histcs1xn`.
//!#     logbase - The logarithm base, if set to `one`, loge is used i.e nat; `two` for bit unit.

	int i       = 0;
	const int c = nbins < 1 ? 1 : nbins;
	double e    = 0.0;

	if (npts > 1 && logbase > 0) {
		for (; i < npts; i++) {
			const double p  = mc_cast(double, h[i]) / c;
			e              -= p > 0.0 ? p * mc_logbase(p, logbase) : 0.0;
		}
	}
	return e;
}

MC_TARGET_FUNC long double mc_histce1xnl(int npts, int nbins, int logbase, const int * h)
{
//!# Requires h[npts] where 1 < n and 0 < nbins and 0 < b. Computing the entropy
//!# of a given `counted` histogram (frequency) @see  `mc_histcs1xn` or `mc_histcg1xn`.
//!#     npts    - Size of histogram i.e h[npts]. @note npts, nbins are usually equal @see `mc_histcs1xn`.
//!#     nbins   - The given bin width.           @note npts, nbins are usually equal @see `mc_histcs1xn`.
//!#     logbase - The logarithm base, if set to `one`, loge is used i.e nat; `two` for bit unit.

	int i         = 0;
	const int c   = nbins < 1 ? 1 : nbins;
	long double e = 0.0L;

	if (npts > 1 && logbase > 0) {
		for (; i < npts; i++) {
			const long double p  = mc_cast(long double, h[i]) / c;
			e                   -= p > 0.0L ? p * mc_logbasel(p, logbase) : 0.0L;
		}
	}
	return e;
}

#endif /* !MC_HISTCE1XN_H */

/* EOF */