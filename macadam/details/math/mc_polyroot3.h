//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_polyroot3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zpolyroot3.h>

#ifndef MC_POLYROOT3_H
#define MC_POLYROOT3_H

#pragma mark - mc_polyroot3 -

MC_TARGET_PROC int mc_polyroot3f(const float a, const float b, const float c, const float d
	, float * r1
	, float * r2
	, float * r3
) {
	int r;
	float i1, i2, i3;
	r = mc_zpolyroot3f(a, b, c, d, r1, &i1, r2, &i2, r3, &i3);
	if (!(r == 2 || r == 3)) {
		*r1 = -1.0f;
		*r2 = -1.0f;
		*r3 = -1.0f;
		r   = -1;
	}
	return r;
}

MC_TARGET_PROC int mc_polyroot3(const double a, const double b, const double c, const double d
	, double * r1
	, double * r2
	, double * r3
) {
	int r;
	double i1, i2, i3;
	r = mc_zpolyroot3(a, b, c, d, r1, &i1, r2, &i2, r3, &i3);
	if (!(r == 2 || r == 3)) {
		*r1 = -1.0;
		*r2 = -1.0;
		*r3 = -1.0;
		r   = -1;
	}
	return r;
}

MC_TARGET_PROC int mc_polyroot3l(const long double a, const long double b, const long double c, const long double d
	, long double * r1
	, long double * r2
	, long double * r3
) {
	int r;
	long double i1, i2, i3;
	r = mc_zpolyroot3l(a, b, c, d, r1, &i1, r2, &i2, r3, &i3);
	if (!(r == 2 || r == 3)) {
		*r1 = -1.0L;
		*r2 = -1.0L;
		*r3 = -1.0L;
		r   = -1;
	}
	return r;
}

#endif /* !MC_POLYROOT3_H */

/* EOF */