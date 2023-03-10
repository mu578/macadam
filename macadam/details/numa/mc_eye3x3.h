//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_eye3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_EYE3X3_H
#define MC_EYE3X3_H

#pragma mark - mc_eye3x3 -

MC_TARGET_FUNC void mc_eye3x3f(float a[9], const int f)
{
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		a[1] = 0.0f; a[2] = 0.0f;
		a[3] = 0.0f; a[5] = 0.0f;
		a[6] = 0.0f; a[7] = 0.0f;
	}
	a[0] = 1.0f; a[4] = 1.0f; a[8] = 1.0f;
}

MC_TARGET_FUNC void mc_eye3x3(double a[9], const int f)
{
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		a[1] = 0.0; a[2] = 0.0;
		a[3] = 0.0; a[5] = 0.0;
		a[6] = 0.0; a[7] = 0.0;
	}
	a[0] = 1.0; a[4] = 1.0; a[8] = 1.0;
}

MC_TARGET_FUNC void mc_eye3x3l(long double a[9], const int f)
{
//!# f=0: set main diagonal to ones and zeroing other elements.
//!# f=1: only set main diagonal to ones.
	if (f == 0) {
		a[1] = 0.0L; a[2] = 0.0L;
		a[3] = 0.0L; a[5] = 0.0L;
		a[6] = 0.0L; a[7] = 0.0L;
	}
	a[0] = 1.0L; a[4] = 1.0L; a[8] = 1.0L;
}

#endif /* !MC_EYE3X3_H */

/* EOF */