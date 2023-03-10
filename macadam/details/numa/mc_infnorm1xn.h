//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_infnorm1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>

#ifndef MC_INFNORM1XN_H
#define MC_INFNORM1XN_H

#pragma mark - mc_infnorm1xn -

MC_TARGET_FUNC float mc_infnorm1xnf(const int n, const float * x, const int f)
{
//!# Requires x[1 x n]. Returning the infinity norm of x.
//!# f=0: computing the maximum of the absolute values.
//!# f=1: computing the minimum of the absolute values.
	int i = 0;
	float nrm = 0.0f, xi;
	for (; i < n ; i++) {
		xi  = mc_fabsf(x[i]);
		nrm = (f == 1) ? ((nrm > xi) ? xi : nrm) : ((nrm > xi) ? nrm : xi);
	}
	return nrm;
}

MC_TARGET_FUNC double mc_infnorm1xnff(const int n, const float * x, const int f)
{
//!# Requires x[1 x n]. Returning the infinity norm of x.
//!# f=0: computing the maximum of the absolute values.
//!# f=1: computing the minimum of the absolute values.
	int i = 0;
	double nrm = 0.0, xi;
	for (; i < n ; i++) {
		xi  = mc_fabs(mc_cast(double, x[i]));
		nrm = (f == 1) ? ((nrm > xi) ? xi : nrm) : ((nrm > xi) ? nrm : xi);
	}
	return nrm;
}

MC_TARGET_FUNC double mc_infnorm1xn(const int n, const double * x, const int f)
{
//!# Requires x[1 x n]. Returning the infinity norm of x.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	int i = 0;
	double nrm = 0.0, xi;
	for (; i < n ; i++) {
		xi  = mc_fabs(x[i]);
		nrm = (f == 1) ? ((nrm > xi) ? xi : nrm) : ((nrm > xi) ? nrm : xi);
	}
	return nrm;
}

MC_TARGET_FUNC long double mc_infnorm1xnl(const int n, const long double * x, const int f)
{
//!# Requires x[1 x n]. Returning the infinity norm of x.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	int i = 0;
	long double nrm = 0.0L, xi;
	for (; i < n ; i++) {
		xi  = mc_fabsl(x[i]);
		nrm = (f == 1) ? ((nrm > xi) ? xi : nrm) : ((nrm > xi) ? nrm : xi);
	}
	return nrm;
}

#endif /* !MC_INFNORM1XN_H */

/* EOF */