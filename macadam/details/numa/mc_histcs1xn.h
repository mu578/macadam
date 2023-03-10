//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_histcs1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/numa/mc_izeros1xn.h>
#include <macadam/details/numa/mc_minmax1xn.h>

#ifndef MC_HISTCS1XN_H
#define MC_HISTCS1XN_H

#pragma mark - mc_histcs1xn -

MC_TARGET_FUNC int mc_histcs1xnf(const int n, const float * x, int adiff, float min, float max, int npts, int nbins, int * h)
{
//!# Requires x[n] and h[npts] where 1 < n.
//!#     n     - Number of samples in x.
//!#     x     - The sample vector `x`.
//!#     adiff - Pass 1 to compute absolute differencials or 0 (default).
//!#     min   - Minimum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     max   - Maximum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     npts  - Size of histogram i.e h[npts].
//!#     nbins - The bin width.
//!#     h     - The histogram result.

	int i = 0;
	float w, scale;

	if (n > 1 && npts > 1) {
		mc_izeros1xn(npts, h);
		if (min == 0.0f && max == 0.0f) {
			mc_minmax1xnf(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		if (nbins < 2) {
			nbins = 2;
		}
		scale = mc_cast_expr(float, (adiff > 0 ? mc_fabsf(max - min) : (max - min)) / (nbins - 1));
		scale = (scale != 0.0f) ? (1.0f / scale) : 1.0f;
		for (; i < n; i++) {
			const int j = mc_cast_expr(const int, (adiff > 0 ? mc_fabsf(x[i] - min) : (x[i] - min)) * scale);
			const int q = j % npts;
 			const int k = q < 0 ? (q + (q >> (MCLIMITS_CBITS * sizeof(int) - 1U) & npts)) : q;
			h[k]++;
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_histcs1xn(const int n, const double * x, int adiff, double min, double max, int npts, int nbins, int * h)
{
//!# Requires x[n] and h[npts] where 1 < n.
//!#     n     - Number of samples in x.
//!#     x     - The sample vector `x`.
//!#     adiff - Pass 1 to compute absolute differencials or 0 (default).
//!#     min   - Minimum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     max   - Maximum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     npts  - Size of histogram i.e h[npts].
//!#     nbins - The bin width.
//!#     h     - The histogram result.

	int i = 0;
	double w, scale;

	if (n > 1 && npts > 1) {
		mc_izeros1xn(npts, h);
		if (min == 0.0 && max == 0.0) {
			mc_minmax1xn(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		if (nbins < 2) {
			nbins = 2;
		}
		scale = mc_cast_expr(double, (adiff > 0 ? mc_fabs(max - min) : (max - min)) / (nbins - 1));
		scale = (scale != 0.0) ? (1.0 / scale) : 1.0;
		for (; i < n; i++) {
			const int j = mc_cast_expr(const int, (adiff > 0 ? mc_fabs(x[i] - min) : (x[i] - min)) * scale);
			const int q = j % npts;
 			const int k = q < 0 ? (q + (q >> (MCLIMITS_CBITS * sizeof(int) - 1U) & npts)) : q;
			h[k]++;
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_histcs1xnl(const int n, const long double * x, int adiff, long double min, long double max, int npts, int nbins, int * h)
{
//!# Requires x[n] and h[npts] where 1 < n.
//!#     n     - Number of samples in x.
//!#     x     - The sample vector `x`.
//!#     adiff - Pass 1 to compute absolute differencials or 0 (default).
//!#     min   - Minimum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     max   - Maximum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     npts  - Size of histogram i.e h[npts].
//!#     nbins - The bin width.
//!#     h     - The histogram result.

	int i = 0;
	long double w, scale;

	if (n > 1 && npts > 1) {
		mc_izeros1xn(npts, h);
		if (min == 0.0L && max == 0.0L) {
			mc_minmax1xnl(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		if (nbins < 2) {
			nbins = 2;
		}
		scale = mc_cast_expr(long double, (adiff > 0 ? mc_fabsl(max - min) : (max - min)) / (nbins - 1));
		scale = (scale != 0.0L) ? (1.0L / scale) : 1.0L;
		for (; i < n; i++) {
			const int j = mc_cast_expr(const int, (adiff > 0 ? mc_fabsl(x[i] - min) : (x[i] - min)) * scale);
			const int q = j % npts;
 			const int k = q < 0 ? (q + (q >> (MCLIMITS_CBITS * sizeof(int) - 1U) & npts)) : q;
			h[k]++;
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_HISTCS1XN_H */

/* EOF */