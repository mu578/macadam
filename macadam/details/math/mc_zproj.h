//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zproj.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_isinf.h>

#ifndef MC_ZPROJ_H
#define MC_ZPROJ_H

#pragma mark - mc_zproj -

MC_TARGET_PROC void mc_zprojf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	*c_r = a_r;
	*c_i = a_i;
	if (mc_isinf(*c_r) || mc_isinf(*c_i)) {
		*c_r = MCK_INFP;
		*c_i = mc_copysignf(0.0f, *c_i);
	}
}

MC_TARGET_PROC void mc_zproj(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	*c_r = a_r;
	*c_i = a_i;
	if (mc_isinf(*c_r) || mc_isinf(*c_i)) {
		*c_r = MCK_INFP;
		*c_i = mc_copysign(0.0, *c_i);
	}
}

MC_TARGET_PROC void mc_zprojl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	*c_r = a_r;
	*c_i = a_i;
	if (mc_isinf(*c_r) || mc_isinf(*c_i)) {
		*c_r = MCK_INFP;
		*c_i = mc_copysignl(0.0L, *c_i);
	}
}

#endif /* !MC_ZPROJ_H */

/* EOF */