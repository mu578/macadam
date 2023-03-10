//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zgamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zexp.h>
#include <macadam/details/math/mc_zgammaln.h>

#ifndef MC_ZGAMMA_H
#define MC_ZGAMMA_H

#pragma mark - mc_zgamma_approx0 -

MC_TARGET_PROC void mc_zgammaf_approx0(float * r_r, float * r_i, float x_r, float x_i)
{
	mc_zgammalnf_approx0(r_r, r_i, x_r, x_i);
	if (*r_r != 0.0f && *r_i != 0.0f) {
		mc_zexpf(r_r, r_i, *r_r, *r_i);
	}
}

MC_TARGET_PROC void mc_zgamma_approx0(double * r_r, double * r_i, double x_r, double x_i)
{
	mc_zgammaln_approx0(r_r, r_i, x_r, x_i);
	if (*r_r != 0.0 && *r_i != 0.0) {
		mc_zexp(r_r, r_i, *r_r, *r_i);
	}
}

MC_TARGET_PROC void mc_zgammal_approx0(long double * r_r, long double * r_i, long double x_r, long double x_i)
{
	mc_zgammalnl_approx0(r_r, r_i, x_r, x_i);
	if (*r_r != 0.0L && *r_i != 0.0L) {
		mc_zexpl(r_r, r_i, *r_r, *r_i);
	}
}

#pragma mark - mc_zgamma -

MC_TARGET_PROC void mc_zgammaf(float * r_r, float * r_i, const float x_r, const float x_i)
{
	mc_zgammaf_approx0(r_r, r_i, x_r, x_i);
}

MC_TARGET_PROC void mc_zgamma(double * r_r, double * r_i, const double x_r, const double x_i)
{
	mc_zgamma_approx0(r_r, r_i, x_r, x_i);
}

MC_TARGET_PROC void mc_zgammal(long double * r_r, long double * r_i, const long double x_r, const long double x_i)
{
	mc_zgammal_approx0(r_r, r_i, x_r, x_i);
}

#endif /* !MC_ZLGAMMA_H */

/* EOF */