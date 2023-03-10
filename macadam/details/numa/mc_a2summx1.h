//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_a2summx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_twosum.h>

#ifndef MC_ASUM2MX1_H
#define MC_ASUM2MX1_H

#pragma mark - mc_a2summx1 -

MC_TARGET_FUNC float mc_a2summx1f(const int m, const int n, const int j, const float * a)
{
	int i   = 1;
	float e = 0.0f, s = 0.0f, y;
	if (n > 0) {
		s = mc_fabsf(a[j]);
		for (; i < m; i++) {
			mc_twosumf(s, mc_fabsf(a[(n * i) + j]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_a2summx1ff(const int m, const int n, const int j, const float * a)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = mc_fabs(mc_cast(double, a[j]));
		for (; i < m; i++) {
			mc_twosum(s, mc_fabs(mc_cast(double, a[(n * i) + j])), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_a2summx1(const int m, const int n, const int j, const double * a)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = mc_fabs(a[j]);
		for (; i < m; i++) {
			mc_twosum(s, mc_fabs(a[(n * i) + j]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC long double mc_a2summx1l(const int m, const int n, const int j, const long double * a)
{
	int i         = 1;
	long double e = 0.0L, s = 0.0L, y;
	if (n > 0) {
		s = mc_fabsl(a[j]);
		for (; i < m; i++) {
			mc_twosuml(s, mc_fabsl(a[(n * i) + j]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

#endif /* !MC_ASUM2MX1_H */

/* EOF */