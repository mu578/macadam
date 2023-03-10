//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sumsq1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/numa/mc_ssqr1xn.h>

#ifndef MC_SUMSQ1XN_H
#define MC_SUMSQ1XN_H

#pragma mark - mc_sumsq1xn -

MC_TARGET_FUNC float mc_sumsq1xnf(const int n, const float * x)
{
	float sumsq = 0.0f, scale;
	if (n > 0) {
		mc_ssqr1xnf(n, x, &sumsq, &scale);
		sumsq = mc_raise2f(scale) * sumsq;
	}
	return sumsq;
}

MC_TARGET_FUNC double mc_sumsq1xnff(const int n, const float * x)
{
	double sumsq = 0.0, scale;
	if (n > 0) {
		mc_ssqr1xnff(n, x, &sumsq, &scale);
		sumsq = mc_raise2(scale) * sumsq;
	}
	return sumsq;
}

MC_TARGET_FUNC double mc_sumsq1xn(const int n, const double * x)
{
	double sumsq = 0.0, scale;
	if (n > 0) {
		mc_ssqr1xn(n, x, &sumsq, &scale);
		sumsq = mc_raise2(scale) * sumsq;
	}
	return sumsq;
}

MC_TARGET_FUNC long double mc_sumsq1xnl(const int n, const long double * x)
{
	long double sumsq = 0.0L, scale;
	if (n > 0) {
		mc_ssqr1xnl(n, x, &sumsq, &scale);
		sumsq = mc_raise2l(scale) * sumsq;
	}
	return sumsq;
}

#endif /* !MC_SUMSQ1XN_H */

/* EOF */