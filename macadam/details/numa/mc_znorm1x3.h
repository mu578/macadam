//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_znorm1x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zmul.h>

#ifndef MC_ZNORM1X3_H
#define MC_ZNORM1X3_H

#pragma mark - mc_znorm1x3 -

MC_TARGET_PROC float mc_znorm1x3f(float x0_r, float x0_i, float x1_r, float x1_i, float x2_r, float x2_i)
{
	float u_r, u_i, w_r, w_i, z_r, z_i;
	mc_zmulf(&u_r, &u_i, x0_r, x0_i, x0_r, -x0_i);
	mc_zmulf(&w_r, &w_i, x1_r, x1_i, x1_r, -x1_i);
	mc_zmulf(&z_r, &z_i, x2_r, x2_i, x2_r, -x2_i);
	return mc_sqrtf(u_r + w_r + z_r);
}

MC_TARGET_PROC double mc_znorm1x3(double x0_r, double x0_i, double x1_r, double x1_i, double x2_r, double x2_i)
{
	double u_r, u_i, w_r, w_i, z_r, z_i;
	mc_zmul(&u_r, &u_i, x0_r, x0_i, x0_r, -x0_i);
	mc_zmul(&w_r, &w_i, x1_r, x1_i, x1_r, -x1_i);
	mc_zmul(&z_r, &z_i, x2_r, x2_i, x2_r, -x2_i);
	return mc_sqrt(u_r + w_r + z_r);
}

MC_TARGET_PROC long double mc_znorm1x3l(long double x0_r, long double x0_i, long double x1_r, long double x1_i, long double x2_r, long double x2_i)
{
	long double u_r, u_i, w_r, w_i, z_r, z_i;
	mc_zmull(&u_r, &u_i, x0_r, x0_i, x0_r, -x0_i);
	mc_zmull(&w_r, &w_i, x1_r, x1_i, x1_r, -x1_i);
	mc_zmull(&z_r, &z_i, x2_r, x2_i, x2_r, -x2_i);
	return mc_sqrtl(u_r + w_r + z_r);
}

#endif /* !MC_ZNORM1X3_H */

/* EOF */