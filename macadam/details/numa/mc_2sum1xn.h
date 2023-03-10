//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_2sum1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_twosum.h>

#ifndef MC_SUM21XN_H
#define MC_SUM21XN_H

#pragma mark - mc_2sum1xn -

MC_TARGET_FUNC float mc_2sum1xnf(const int n, const float * x)
{
	int i   = 1;
	float e = 0.0f, s = 0.0f, y;
	if (n > 0) {
		s = x[0];
		for (; i < n; i++) {
			mc_twosumf(s, x[i], &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_2sum1xnff(const int n, const float * x)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = mc_cast(double, x[0]);
		for (; i < n; i++) {
			mc_twosum(s, mc_cast(double, x[i]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_2sum1xn(const int n, const double * x)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = x[0];
		for (; i < n; i++) {
			mc_twosum(s, x[i], &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC long double mc_2sum1xnl(const int n, const long double * x)
{
	int i         = 1;
	long double e = 0.0L, s = 0.0L, y;
	if (n > 0) {
		s = x[0];
		for (; i < n; i++) {
			mc_twosuml(s, x[i], &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

#endif /* !MC_SUM21XN_H */

/* EOF */