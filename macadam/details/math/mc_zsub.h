//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zsub.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ZSUB_H
#define MC_ZSUB_H

#pragma mark - mc_zsub -

MC_TARGET_PROC void mc_zsubf(float * c_r, float * c_i
	, const float a_r, const float a_i
	, const float b_r, const float b_i
) {
	*c_r = a_r - b_r;
	*c_i = a_i - b_i;
}

MC_TARGET_PROC void mc_zsub(double * c_r, double * c_i
	, const double a_r, const double a_i
	, const double b_r, const double b_i
) {
	*c_r = a_r - b_r;
	*c_i = a_i - b_i;
}

MC_TARGET_PROC void mc_zsubl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
	, const long double b_r, const long double b_i
) {
	*c_r = a_r - b_r;
	*c_i = a_i - b_i;
}

#endif /* !MC_ZSUB_H */

/* EOF */