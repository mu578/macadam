//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_polyroot2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zpolyroot2.h>

#ifndef MC_POLYROOT2_H
#define MC_POLYROOT2_H

#pragma mark - mc_polyroot2 -

MC_TARGET_PROC int mc_polyroot2f(const float a, const float b, const float c
	, float * r1
	, float * r2
) {
	int r;
	float i1, i2;
	r = mc_zpolyroot2f(a, b, c, r1, &i1, r2, &i2);
	if (!(r == 1 || r == 2)) {
		*r1 = -1.0f;
		*r2 = -1.0f;
		r   = -1;
	}
	return r;
}

MC_TARGET_PROC int mc_polyroot2(const double a, const double b, const double c
	, double * r1
	, double * r2
) {
	int r;
	double i1, i2;
	r = mc_zpolyroot2(a, b, c, r1, &i1, r2, &i2);
	if (!(r == 1 || r == 2)) {
		*r1 = -1.0;
		*r2 = -1.0;
		r   = -1;
	}
	return r;
}

MC_TARGET_PROC int mc_polyroot2l(const long double a, const long double b, const long double c
	, long double * r1
	, long double * r2
) {
	int r;
	long double i1, i2;
	r = mc_zpolyroot2l(a, b, c, r1, &i1, r2, &i2);
	if (!(r == 1 || r == 2)) {
		*r1 = -1.0L;
		*r2 = -1.0L;
		r   = -1;
	}
	return r;
}

#endif /* !MC_POLYROOT2_H */

/* EOF */