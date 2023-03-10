//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_histcg1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_izeros1xn.h>

#ifndef MC_HISTCG1XN_H
#define MC_HISTCG1XN_H

#pragma mark - mc_histcg1xn -

MC_TARGET_FUNC int mc_histcg1xnf(const int n, const float * x, float min, float max, int nbins, int * h)
{
//!# Requires x[n] and h[nbins] where 1 < n && 0 < nbins.
//!#     n     - Number of samples in x.
//!#     x     - The sample vector `x`.
//!#     min   - Minimum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     max   - Maximum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     nbins - The given bin width.
//!#     h     - The histogram result.

	int i = 0;
	float w;

	if (n > 1 && nbins > 0) {
		if (min == 0.0f && max == 0.0f) {
			mc_minmax1xnf(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		if (nbins < 0x1000001) {
			const float c = mc_cast_expr(const float, (max - min) / nbins);
			if (c != 0.0f) {
				mc_izeros1xn(nbins, h);
				for (; i < n; i++) {
					const int j = mc_cast_expr(const int, (x[i] - min) / c);
					if (j >= 0 && j < nbins) {
						h[j]++;
					}
				}
				return 0;
			}
		} else {
			const double c = mc_cast_expr(const double, (max - min) / nbins);
			if (c != 0.0) {
				mc_izeros1xn(nbins, h);
				for (; i < n; i++) {
					const int j = mc_cast_expr(const int, mc_cast_expr(const double, (x[i] - min)) / c);
					if (j >= 0 && j < nbins) {
						h[j]++;
					}
				}
				return 0;
			}
		}
	}
	return -1;
}

MC_TARGET_FUNC int mc_histcg1xn(const int n, const double * x, double min, double max, int nbins, int * h)
{
//!# Requires x[n] and h[nbins] where 1 < n && 0 < nbins.
//!#     n     - Number of samples in x.
//!#     x     - The sample vector `x`.
//!#     min   - Minimum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     max   - Maximum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     nbins - The given bin width.
//!#     h     - The histogram result.

	int i = 0;
	double w;

	if (n > 1 && nbins > 0) {
		if (min == 0.0f && max == 0.0f) {
			mc_minmax1xn(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		const double c = mc_cast_expr(const double, (max - min) / nbins);
		if (c != 0.0) {
			mc_izeros1xn(nbins, h);
			for (; i < n; i++) {
				const int j = mc_cast_expr(const int, (x[i] - min) / c);
				if (j >= 0 && j < nbins) {
					h[j]++;
				}
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC  int mc_histcg1xnl(const int n, const long double * x, long double min, long double max, int nbins, int * h)
{
//!# Requires x[n] and h[nbins] where 1 < n && 0 < nbins.
//!#     n     - Number of samples in x.
//!#     x     - The sample vector `x`.
//!#     nbins - The given bin width.
//!#     min   - Minimum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     max   - Maximum edge value. If min and max are set to `zero` computing true x-min and x-max.
//!#     h     - The histogram result.

	int i = 0;
	long double w;

	if (n > 1 && nbins > 0) {
		if (min == 0.0f && max == 0.0f) {
			mc_minmax1xnl(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		const long double c = mc_cast_expr(const long double, (max - min) / nbins);
		if (c != 0.0) {
			mc_izeros1xn(nbins, h);
			for (; i < n; i++) {
				const int j = mc_cast_expr(const int, (x[i] - min) / c);
				if (j >= 0 && j < nbins) {
					h[j]++;
				}
			}
			return 0;
		}
	}
	return -1;
}

#endif /* !MC_HISTCG1XN_H */

/* EOF */