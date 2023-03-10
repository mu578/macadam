//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zeye2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ZEYE2X2_H
#define MC_ZEYE2X2_H

#pragma mark - mc_zeye2x2 -

MC_TARGET_PROC void mc_zeye2x2f(
	  float * v0_r, float * v0_i, float * v1_r, float * v1_i
	, float * v2_r, float * v2_i, float * v3_r, float * v3_i
	, const int f
) {
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		*v1_r = 0.0f; *v1_i = 0.0f;
		*v2_r = 0.0f; *v2_i = 0.0f;
	}
	*v0_r = 1.0f; *v0_i = 0.0f; *v3_r = 1.0f; *v3_i = 0.0f;
}

MC_TARGET_PROC void mc_zeye2x2(
	  double * v0_r, double * v0_i, double * v1_r, double * v1_i
	, double * v2_r, double * v2_i, double * v3_r, double * v3_i
	, const int f
) {
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		*v1_r = 0.0; *v1_i = 0.0;
		*v2_r = 0.0; *v2_i = 0.0;
	}
	*v0_r = 1.0; *v0_i = 0.0; *v3_r = 1.0; *v3_i = 0.0;
}

MC_TARGET_PROC void mc_zeye2x2l(
	  long double * v0_r, long double * v0_i, long double * v1_r, long double * v1_i
	, long double * v2_r, long double * v2_i, long double * v3_r, long double * v3_i
	, const int f
) {
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		*v1_r = 0.0L; *v1_i = 0.0L;
		*v2_r = 0.0L; *v2_i = 0.0L;
	}
	*v0_r = 1.0L; *v0_i = 0.0L; *v3_r = 1.0L; *v3_i = 0.0L;
}

#endif /* !MC_ZEYE2X2_H */

/* EOF */