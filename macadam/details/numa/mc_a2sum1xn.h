//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_a2sum1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_twosum.h>

#ifndef MC_ASUM21XN_H
#define MC_ASUM21XN_H

#pragma mark - mc_a2sum1xn -

MC_TARGET_FUNC float mc_a2sum1xnf(const int n, const float * x)
{
	int i   = 1;
	float e = 0.0f, s = 0.0f, y;
	if (n > 0) {
		s = mc_fabsf(x[0]);
		for (; i < n; i++) {
			mc_twosumf(s, mc_fabsf(x[i]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_a2sum1xnff(const int n, const float * x)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = mc_fabs(mc_cast(double, x[0]));
		for (; i < n; i++) {
			mc_twosum(s, mc_fabs(mc_cast(double, x[i])), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_a2sum1xn(const int n, const double * x)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = mc_fabs(x[0]);
		for (; i < n; i++) {
			mc_twosum(s, mc_fabs(x[i]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC long double mc_a2sum1xnl(const int n, const long double * x)
{
	int i         = 1;
	long double e = 0.0L, s = 0.0L, y;
	if (n > 0) {
		s = mc_fabsl(x[0]);
		for (; i < n; i++) {
			mc_twosuml(s, mc_fabsl(x[i]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

#endif /* !MC_ASUM21XN_H */

/* EOF */