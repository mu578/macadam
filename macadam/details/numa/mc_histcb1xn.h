//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_histcb1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cbrti.h>
#include <macadam/details/math/mc_ceil.h>
#include <macadam/details/numa/mc_minmax1xn.h>
#include <macadam/details/numa/mc_stdd1xn.h>

#ifndef MC_HISTNB1XN_H
#define MC_HISTNB1XN_H

#pragma mark - mc_histcb1xn -

MC_TARGET_FUNC int mc_histcb1xnf(const int n, const float * x, float min, float max, float stddev, float tol)
{
	float w;

	if (n > 1) {
		if (min == 0.0f && max == 0.0f) {
			mc_minmax1xnf(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		if (stddev == 0.0f) {
			mc_stdd1xnf(n, x, 1, &stddev);
		}
		if (tol <= 0.0f) {
			tol = 1E-04f;
		}
		if (stddev > tol) {
			w = mc_ceilf((min - max) / (3.5f * stddev / mc_cbrtif(n)));
			return mc_cast(int, w);
		}
	}
	return -1;
}

MC_TARGET_FUNC int mc_histcb1xn(const int n, const double * x, double min, double max, double stddev, double tol)
{
	double w;

	if (n > 1) {
		if (min == 0.0 && max == 0.0) {
			mc_minmax1xn(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		if (stddev == 0.0) {
			mc_stdd1xn(n, x, 1, &stddev);
		}
		if (tol <= 0.0) {
			tol = 1E-04;
		}
		if (stddev > tol) {
			w = mc_ceil((min - max) / (3.5 * stddev / mc_cbrti(n)));
			return mc_cast(int, w);
		}
	}
	return -1;
}

MC_TARGET_FUNC int mc_histcb1xnl(const int n, const long double * x, long double min, long double max, long double stddev, long double tol)
{
	long double w;

	if (n > 1) {
		if (min == 0.0L && max == 0.0L) {
			mc_minmax1xnl(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
		} else if (min > max) {
			mcswap_var(w, min, max);
		}
		if (stddev == 0.0L) {
			mc_stdd1xnl(n, x, 1, &stddev);
		}
		if (tol <= 0.0L) {
			tol = 1E-04L;
		}
		if (stddev > tol) {
			w = mc_ceill((min - max) / (3.5L * stddev / mc_cbrtil(n)));
			return mc_cast(int, w);
		}
	}
	return -1;
}

#endif /* !MC_HISTNB1XN_H */

/* EOF */