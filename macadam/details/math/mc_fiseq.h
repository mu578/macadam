//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fiseq.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_absmag.h>
#include <macadam/details/math/mc_fabs.h>

#ifndef MC_FISEQ_H
#define MC_FISEQ_H

#pragma mark - mc_fiseq -

MC_TARGET_FUNC int mc_fiseqf(const float a, const float b, const float p, const int n)
{
//!# Performing real-floating-point equivalence between `a` and `b` at `n` arbitrary epsilon
//!# precision denoted `p`. If `p` is set to `zero`, `p` is reset to the default machine-epsilon.
	float d, e;

	if (a == b) {
		return 1;
	} else {
		const float eps = (
			  mc_iabs(n) != 0 || mc_iabs(n) != 1
			? mc_cast(float, mc_iabs(n)) * (p == 0.0f ? MCLIMITS_EPSILONF : p)
			: (p == 0.0f ? MCLIMITS_EPSILONF : p)
		);
		if (a == 0.0f) {
			if (mc_fabsf(b) <= eps) {
				return 1;
			}
		} else if (b == 0.0f) {
			if (mc_fabsf(a) <= eps) {
				return 1;
			}
		} else {
			e = 2.0f * a * eps;
			d = a - b;
			if (e < 0.0f) {
				e = -e;
			}
			if (-e <= d && d <= e) {
				return 1;
			}
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_fiseq(const double a, const double b, const float p, const int n)
{
//!# Performing real-floating-point equivalence between `a` and `b` at `n` arbitrary epsilon
//!# precision denoted `p`. If `p` is set to `zero`, `p` is reset to the default machine-epsilon.
	double d, e;

	if (a == b) {
		return 1;
	} else {
		const double eps = (
			  mc_iabs(n) != 0 || mc_iabs(n) != 1
			? mc_cast(double, mc_iabs(n)) * (p == 0.0 ? MCLIMITS_EPSILON : p)
			: (p == 0.0 ? MCLIMITS_EPSILON : p)
		);
		if (a == 0.0) {
			if (mc_fabs(b) <= eps) {
				return 1;
			}
		} else if (b == 0.0) {
			if (mc_fabs(a) <= eps) {
				return 1;
			}
		} else {
			e = 2.0 * a * eps;
			d = a - b;
			if (e < 0.0) {
				e = -e;
			}
			if (-e <= d && d <= e) {
				return 1;
			}
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_fiseql(const long double a, const long double b, const float p, const int n)
{
//!# Performing real-floating-point equivalence between `a` and `b` at `n` arbitrary epsilon
//!# precision denoted `p`. If `p` is set to `zero`, `p` is reset to the default machine-epsilon.
	long double d, e;

	if (a == b) {
		return 1;
	} else {
		const long double eps = (
			  mc_iabs(n) != 0 || mc_iabs(n) != 1
			? mc_cast(long double, mc_iabs(n)) * (p == 0.0L ? MCLIMITS_EPSILONL : p)
			: (p == 0.0L ? MCLIMITS_EPSILONL : p)
		);
		if (a == 0.0L) {
			if (mc_fabsl(b) <= eps) {
				return 1;
			}
		} else if (b == 0.0L) {
			if (mc_fabsl(a) <= eps) {
				return 1;
			}
		} else {
			e = 2.0L * a * eps;
			d = a - b;
			if (e < 0.0L) {
				e = -e;
			}
			if (-e <= d && d <= e) {
				return 1;
			}
		}
	}
	return 0;
}

#endif /* !MC_FISEQ_H */

/* EOF */