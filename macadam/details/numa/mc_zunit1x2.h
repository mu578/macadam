//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zunit1x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zdiv.h>
#include <macadam/details/numa/mc_znorm1x2.h>

#ifndef MC_ZUNIT1X2_H
#define MC_ZUNIT1X2_H

#pragma mark - mc_zunit1x2 -

MC_TARGET_FUNC void mc_zunit1x2f(float * x0_r, float * x0_i, float * x1_r, float * x1_i)
{
	const float nrm2 = mc_znorm1x2f(*x0_r, *x0_i, *x1_r, *x1_i);
	mc_zdivf(x0_r, x0_i, *x0_r, *x0_i, nrm2, 0.0f);
	mc_zdivf(x1_r, x1_i, *x1_r, *x1_i, nrm2, 0.0f);
}

MC_TARGET_FUNC void mc_zunit1x2(double * x0_r, double * x0_i, double * x1_r, double * x1_i)
{
	const double nrm2 = mc_znorm1x2(*x0_r, *x0_i, *x1_r, *x1_i);
	mc_zdiv(x0_r, x0_i, *x0_r, *x0_i, nrm2, 0.0);
	mc_zdiv(x1_r, x1_i, *x1_r, *x1_i, nrm2, 0.0);
}

MC_TARGET_FUNC void mc_zunit1x2l(long double * x0_r, long double * x0_i, long double * x1_r, long double * x1_i)
{
	const long double nrm2 = mc_znorm1x2l(*x0_r, *x0_i, *x1_r, *x1_i);
	mc_zdivl(x0_r, x0_i, *x0_r, *x0_i, nrm2, 0.0L);
	mc_zdivl(x1_r, x1_i, *x1_r, *x1_i, nrm2, 0.0L);
}

#endif /* !MC_ZUNIT1X2_H */

/* EOF */