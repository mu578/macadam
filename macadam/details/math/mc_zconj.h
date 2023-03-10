//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zconj.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ZCONJ_H
#define MC_ZCONJ_H

#pragma mark - mc_zconj -

MC_TARGET_PROC void mc_zconjf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	*c_r =  a_r;
	*c_i = -a_i;
}

MC_TARGET_PROC void mc_zconj(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	*c_r =  a_r;
	*c_i = -a_i;
}

MC_TARGET_PROC void mc_zconjl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	*c_r =  a_r;
	*c_i = -a_i;
}

#endif /* !MC_ZCONJ_H */

/* EOF */