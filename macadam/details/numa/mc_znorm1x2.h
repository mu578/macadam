//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_znorm1x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zmul.h>

#ifndef MC_ZNORM1X2_H
#define MC_ZNORM1X2_H

#pragma mark - mc_znorm1x2 -

MC_TARGET_PROC float mc_znorm1x2f(float x0_r, float x0_i, float x1_r, float x1_i)
{
	float w_r, w_i, z_r, z_i;
	mc_zmulf(&w_r, &w_i, x0_r, x0_i, x0_r, -x0_i);
	mc_zmulf(&z_r, &z_i, x1_r, x1_i, x1_r, -x1_i);
	return mc_sqrtf(w_r + z_r);
}

MC_TARGET_PROC double mc_znorm1x2(double x0_r, double x0_i, double x1_r, double x1_i)
{
	double w_r, w_i, z_r, z_i;
	mc_zmul(&w_r, &w_i, x0_r, x0_i, x0_r, -x0_i);
	mc_zmul(&z_r, &z_i, x1_r, x1_i, x1_r, -x1_i);
	return mc_sqrt(w_r + z_r);
}

MC_TARGET_PROC long double mc_znorm1x2l(long double x0_r, long double x0_i, long double x1_r, long double x1_i)
{
	long double w_r, w_i, z_r, z_i;
	mc_zmull(&w_r, &w_i, x0_r, x0_i, x0_r, -x0_i);
	mc_zmull(&z_r, &z_i, x1_r, x1_i, x1_r, -x1_i);
	return mc_sqrtl(w_r + z_r);
}

#endif /* !MC_ZNORM1X2_H */

/* EOF */