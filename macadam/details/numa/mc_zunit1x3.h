//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zunit1x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zdiv.h>
#include <macadam/details/numa/mc_znorm1x3.h>

#ifndef MC_ZUNIT1X3_H
#define MC_ZUNIT1X3_H

#pragma mark - mc_zunit1x3 -

MC_TARGET_PROC void mc_zunit1x3f(float * x0_r, float * x0_i, float * x1_r, float * x1_i, float * x2_r, float * x2_i)
{
	const float nrm2 = mc_znorm1x3f(*x0_r, *x0_i, *x1_r, *x1_i, *x2_r, *x2_i);
	mc_zdivf(x0_r, x0_i, *x0_r, *x0_i, nrm2, 0.0f);
	mc_zdivf(x1_r, x1_i, *x1_r, *x1_i, nrm2, 0.0f);
	mc_zdivf(x2_r, x2_i, *x2_r, *x2_i, nrm2, 0.0f);
}

MC_TARGET_PROC void mc_zunit1x3(double * x0_r, double * x0_i, double * x1_r, double * x1_i, double * x2_r, double * x2_i)
{
	const double nrm2 = mc_znorm1x3(*x0_r, *x0_i, *x1_r, *x1_i, *x2_r, *x2_i);
	mc_zdiv(x0_r, x0_i, *x0_r, *x0_i, nrm2, 0.0);
	mc_zdiv(x1_r, x1_i, *x1_r, *x1_i, nrm2, 0.0);
	mc_zdiv(x2_r, x2_i, *x2_r, *x2_i, nrm2, 0.0);
}

MC_TARGET_PROC void mc_zunit1x3l(long double * x0_r, long double * x0_i, long double * x1_r, long double * x1_i, long double * x2_r, long double * x2_i)
{
	const long double nrm2 = mc_znorm1x3l(*x0_r, *x0_i, *x1_r, *x1_i, *x2_r, *x2_i);
	mc_zdivl(x0_r, x0_i, *x0_r, *x0_i, nrm2, 0.0L);
	mc_zdivl(x1_r, x1_i, *x1_r, *x1_i, nrm2, 0.0L);
	mc_zdivl(x2_r, x2_i, *x2_r, *x2_i, nrm2, 0.0L);
}

#endif /* !MC_ZUNIT1X3_H */

/* EOF */