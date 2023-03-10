//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_l2norm1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_ssqr1xn.h>

#ifndef MC_L2NORM1XN_H
#define MC_L2NORM1XN_H

#pragma mark - mc_l2norm1xn -

MC_TARGET_FUNC float mc_l2norm1xnf(const int n, const float * x)
{
	float sumsq = 0.0f, scale;
	if (n > 0) {
		mc_ssqr1xnf(n, x, &sumsq, &scale);
		sumsq = scale * mc_sqrtf(sumsq);
	}
	return sumsq;
}

MC_TARGET_FUNC double mc_l2norm1xnff(const int n, const float * x)
{
	double sumsq = 0.0, scale;
	if (n > 0) {
		mc_ssqr1xnff(n, x, &sumsq, &scale);
		sumsq = scale * mc_sqrt(sumsq);
	}
	return sumsq;
}

MC_TARGET_FUNC double mc_l2norm1xn(const int n, const double * x)
{
	double sumsq = 0.0, scale;
	if (n > 0) {
		mc_ssqr1xn(n, x, &sumsq, &scale);
		sumsq = scale * mc_sqrt(sumsq);
	}
	return sumsq;
}

MC_TARGET_FUNC long double mc_l2norm1xnl(const int n, const long double * x)
{
	long double sumsq = 0.0L, scale;
	if (n > 0) {
		mc_ssqr1xnl(n, x, &sumsq, &scale);
		sumsq = scale * mc_sqrtl(sumsq);
	}
	return sumsq;
}

#endif /* !MC_L2NORM1XN_H */

/* EOF */