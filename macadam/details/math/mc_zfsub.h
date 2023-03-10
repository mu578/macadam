//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zfsub.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ZFSUB_H
#define MC_ZFSUB_H

#pragma mark - mc_zfsub -

MC_TARGET_PROC void mc_zfsubf(float * c_r, float * c_i
	, const float a_r, const float a_i
	, const float b
) {
	*c_r = a_r - b;
	*c_i = a_i;
}

MC_TARGET_PROC void mc_zfsub(double * c_r, double * c_i
	, const double a_r, const double a_i
	, const double b
) {
	*c_r = a_r - b;
	*c_i = a_i;
}

MC_TARGET_PROC void mc_zfsubl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
	, const long double b
) {
	*c_r = a_r - b;
	*c_i = a_i;
}

#endif /* !MC_ZFSUB_H */

/* EOF */