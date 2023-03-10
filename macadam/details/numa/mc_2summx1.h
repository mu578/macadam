//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_2summx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_twosum.h>

#ifndef MC_SUM2MX1_H
#define MC_SUM2MX1_H

#pragma mark - mc_2summx1 -

MC_TARGET_FUNC float mc_2summx1f(const int m, const int n, const int j, const float * a)
{
	int i   = 1;
	float e = 0.0f, s = 0.0f, y;
	if (n > 0) {
		s = a[j];
		for (; i < m; i++) {
			mc_twosumf(s, a[(n * i) + j], &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_2summx1ff(const int m, const int n, const int j, const float * a)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = mc_cast(double, a[j]);
		for (; i < m; i++) {
			mc_twosum(s, mc_cast(double, a[(n * i) + j]), &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC double mc_2summx1(const int m, const int n, const int j, const double * a)
{
	int i    = 1;
	double e = 0.0, s = 0.0, y;
	if (n > 0) {
		s = a[j];
		for (; i < m; i++) {
			mc_twosum(s, a[(n * i) + j], &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

MC_TARGET_FUNC long double mc_2summx1l(const int m, const int n, const int j, const long double * a)
{
	int i         = 1;
	long double e = 0.0L, s = 0.0L, y;
	if (n > 0) {
		s = a[j];
		for (; i < m; i++) {
			mc_twosuml(s, a[(n * i) + j], &s, &y);
			e = e + y;
		}
	}
	return s + e;
}

#endif /* !MC_SUM2MX1_H */

/* EOF */