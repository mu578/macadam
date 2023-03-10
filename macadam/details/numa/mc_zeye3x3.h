//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zeye3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ZEYE3X3_H
#define MC_ZEYE3X3_H

#pragma mark - mc_zeye3x3 -

MC_TARGET_PROC void mc_zeye3x3f(
	  float * v0_r, float * v0_i, float * v1_r, float * v1_i, float * v2_r, float * v2_i
	, float * v3_r, float * v3_i, float * v4_r, float * v4_i, float * v5_r, float * v5_i
	, float * v6_r, float * v6_i, float * v7_r, float * v7_i, float * v8_r, float * v8_i
	, const int f
) {
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		*v1_r = 0.0f; *v1_i = 0.0f; *v2_r = 0.0f; *v2_i = 0.0f;
		*v3_r = 0.0f; *v3_i = 0.0f; *v5_r = 0.0f; *v5_i = 0.0f;
		*v6_r = 0.0f; *v6_i = 0.0f; *v7_r = 0.0f; *v7_i = 0.0f;
	}
	*v0_r = 1.0f; *v0_i = 0.0f; *v4_r = 1.0f; *v4_i = 0.0f; *v8_r = 1.0f; *v8_i = 0.0f;
}

MC_TARGET_PROC void mc_zeye3x3(
	  double * v0_r, double * v0_i, double * v1_r, double * v1_i, double * v2_r, double * v2_i
	, double * v3_r, double * v3_i, double * v4_r, double * v4_i, double * v5_r, double * v5_i
	, double * v6_r, double * v6_i, double * v7_r, double * v7_i, double * v8_r, double * v8_i
	, const int f
) {
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		*v1_r = 0.0; *v1_i = 0.0; *v2_r = 0.0; *v2_i = 0.0;
		*v3_r = 0.0; *v3_i = 0.0; *v5_r = 0.0; *v5_i = 0.0;
		*v6_r = 0.0; *v6_i = 0.0; *v7_r = 0.0; *v7_i = 0.0;
	}
	*v0_r = 1.0; *v0_i = 0.0; *v4_r = 1.0; *v4_i = 0.0; *v8_r = 1.0; *v8_i = 0.0;
}

MC_TARGET_PROC void mc_zeye3x3l(
	  long double * v0_r, long double * v0_i, long double * v1_r, long double * v1_i, long double * v2_r, long double * v2_i
	, long double * v3_r, long double * v3_i, long double * v4_r, long double * v4_i, long double * v5_r, long double * v5_i
	, long double * v6_r, long double * v6_i, long double * v7_r, long double * v7_i, long double * v8_r, long double * v8_i
	, const int f
) {
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		*v1_r = 0.0L; *v1_i = 0.0L; *v2_r = 0.0L; *v2_i = 0.0L;
		*v3_r = 0.0L; *v3_i = 0.0L; *v5_r = 0.0L; *v5_i = 0.0L;
		*v6_r = 0.0L; *v6_i = 0.0L; *v7_r = 0.0L; *v7_i = 0.0L;
	}
	*v0_r = 1.0L; *v0_i = 0.0L; *v4_r = 1.0L; *v4_i = 0.0L; *v8_r = 1.0L; *v8_i = 0.0L;
}

#endif /* !MC_ZEYE3X3_H */

/* EOF */