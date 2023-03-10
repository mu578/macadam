//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mean1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_sum1xn.h>

#ifndef MC_MEAN1XN_H
#define MC_MEAN1XN_H

#pragma mark - mc_mean1xn -

MC_TARGET_FUNC float mc_mean1xnf(const int n, const float * x, const int b, const int f)
{
	float s = 0.0f;
	if (n > 0) {
		s = mc_sum1xnf(n, x, f);
		s = s / mc_cast(const float, (b ? n - 1 : n));
	}
	return s;
}

MC_TARGET_FUNC double mc_mean1xnff(const int n, const float * x, const int b, const int f)
{
	double s = 0.0;
	if (n > 0) {
		s = mc_sum1xnff(n, x, f);
		s = s / mc_cast(const double, (b ? n - 1 : n));
	}
	return s;
}

MC_TARGET_FUNC double mc_mean1xn(const int n, const double * x, const int b, const int f)
{
	double s = 0.0;
	if (n > 0) {
		s = mc_sum1xn(n, x, f);
		s = s / mc_cast(const double, (b ? n - 1 : n));
	}
	return s;
}

MC_TARGET_FUNC long double mc_mean1xnl(const int n, const long double * x, const int b, const int f)
{
	long double s = 0.0L;
	if (n > 0) {
		s = mc_sum1xnl(n, x, f);
		s = s / mc_cast(const long double, (b ? n - 1 : n));
	}
	return s;
}

#endif /* !MC_MEAN1XN_H */

/* EOF */