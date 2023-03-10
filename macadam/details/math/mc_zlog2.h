//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zlog2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zfmul.h>
#include <macadam/details/math/mc_zlog.h>

#ifndef MC_ZLOG2_H
#define MC_ZLOG2_H

#pragma mark - mc_zlog2 -

MC_TARGET_PROC void mc_zlog2f(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	mc_zlogf(c_r, c_i, a_r, a_i);
	mc_zfmulf(c_r, c_i, *c_r, *c_i, MCK_KF(MCK_1_LOGE2));
}

MC_TARGET_PROC void mc_zlog2(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_zlog(c_r, c_i, a_r, a_i);
	mc_zfmul(c_r, c_i, *c_r, *c_i, MCK_K(MCK_1_LOGE2));
}

MC_TARGET_PROC void mc_zlog2l(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_zlogl(c_r, c_i, a_r, a_i);
	mc_zfmull(c_r, c_i, *c_r, *c_i, MCK_KL(MCK_1_LOGE2));
}

#endif /* !MC_ZLOG2_H */

/* EOF */