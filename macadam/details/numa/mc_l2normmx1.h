//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_l2normmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_ssqrmx1.h>

#ifndef MC_L2NORMMX1_H
#define MC_L2NORMMX1_H

#pragma mark - mc_l2normmx1 -

MC_TARGET_FUNC float mc_l2normmx1f(const int m, const int n, const int j, const float * a)
{
	float sumsq = 0.0f, scale;
	if (n > 0) {
		mc_ssqrmx1f(m, n, j, a, &sumsq, &scale);
		sumsq = scale * mc_sqrtf(sumsq);
	}
	return sumsq;
}

MC_TARGET_FUNC double mc_l2normmx1ff(const int m, const int n, const int j, const float * a)
{
	double sumsq = 0.0, scale;
	if (n > 0) {
		mc_ssqrmx1ff(m, n, j, a, &sumsq, &scale);
		sumsq = scale * mc_sqrt(sumsq);
	}
	return sumsq;
}

MC_TARGET_FUNC double mc_l2normmx1(const int m, const int n, const int j, const double * a)
{
	double sumsq = 0.0, scale;
	if (n > 0) {
		mc_ssqrmx1(m, n, j, a, &sumsq, &scale);
		sumsq = scale * mc_sqrt(sumsq);
	}
	return sumsq;
}

MC_TARGET_FUNC long double mc_l2normmx1l(const int m, const int n, const int j, const long double * a)
{
	long double sumsq = 0.0L, scale;
	if (n > 0) {
		mc_ssqrmx1l(m, n, j, a, &sumsq, &scale);
		sumsq = scale * mc_sqrtl(sumsq);
	}
	return sumsq;
}

#endif /* !MC_L2NORMMX1_H */

/* EOF */