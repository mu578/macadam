//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zpolyroot2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zadd.h>
#include <macadam/details/math/mc_zfdiv.h>
#include <macadam/details/math/mc_zsqrt.h>
#include <macadam/details/math/mc_zsub.h>

#ifndef MC_ZPOLYROOT2_H
#define MC_ZPOLYROOT2_H

#pragma mark - mc_zpolyroot2 -

MC_TARGET_PROC int mc_zpolyroot2f(const float a, const float b, const float c
	, float * r1
	, float * i1
	, float * r2
	, float * i2
) {
//!# Two roots of the quadratic equation: ax^2+bx+c=0.
//!# Returns 1 if roots are real and equal.
//!# Returns 2 if roots are real and distinct.
//!# Returns 4 if roots are conjugate.
	int r         = -1;
	const float e = 1E-05f;
	const float s = mc_raise2f(b) - 4.0f * a * c;

	if (mc_fabsf(s) < e) {
		*r1 = -0.5f * b / a;
		*i1 = 0.0f;
		*r2 = *r1;
		*i2 = *i1;
		 r  = 1;
	} else if (s > 0.0f) {
		const float q = mc_sqrtf(s);
		*r1           = (-b + q) * (1.0f / (a * 2.0f));
		*i1           = 0.0f;
		*r2           = (-b - q) * (1.0f / (a * 2.0f));
		*i2           = 0.0f;
		 r            = 2;
	} else if (s < 0.0f) {
		mc_zsqrtf(r1, i1, s, 0.0f);
		*r2 = *r1;
		*i2 = *i1;
		 r  = 4;

		mc_zaddf(r1, i1, -b, 0.0f, *r1, *i1);
		mc_zfdivf(r1, i1, *r1, *i1, a * 2.0f);

		mc_zsubf(r2, i2, -b, 0.0f, *r2, *i2);
		mc_zfdivf(r2, i2, *r2, *i2, a * 2.0f);
	}
	return r;
}

MC_TARGET_PROC int mc_zpolyroot2(double a, double b, double c
	, double * r1
	, double * i1
	, double * r2
	, double * i2
) {
//!# Two roots of the quadratic equation: ax^2+bx+c=0.
//!# Returns 1 if roots are real and equal.
//!# Returns 2 if roots are real and distinct.
//!# Returns 4 if roots are conjugate.
	int r          = -1;
	const double e = 1E-09;
	const double s = mc_raise2(b) - 4.0 * a * c;

	if (mc_fabs(s) < e) {
		*r1 = -0.5 * b / a;
		*i1 = 0.0;
		*r2 = *r1;
		*i2 = *i1;
		 r  = 1;
	} else if (s > 0.0) {
		const double q = mc_sqrt(s);
		*r1            = (-b + q) * (1.0 / (a * 2.0));
		*i1            = 0.0;
		*r2            = (-b - q) * (1.0 / (a * 2.0));
		*i2            = 0.0;
		 r             = 2;
	} else if (s < 0.0) {
		mc_zsqrt(r1, i1, s, 0.0);
		*r2 = *r1;
		*i2 = *i1;
		 r  = 4;

		mc_zadd(r1, i1, -b, 0.0, *r1, *i1);
		mc_zfdiv(r1, i1, *r1, *i1, a * 2.0);

		mc_zsub(r2, i2, -b, 0.0, *r2, *i2);
		mc_zfdiv(r2, i2, *r2, *i2, a * 2.0);
	}
	return r;
}

MC_TARGET_PROC int mc_zpolyroot2l(long double a, long double b, long double c
	, long double * r1
	, long double * i1
	, long double * r2
	, long double * i2
) {
//!# Two roots of the quadratic equation: ax^2+bx+c=0.
//!# Returns 1 if roots are real and equal.
//!# Returns 2 if roots are real and distinct.
//!# Returns 4 if roots are conjugate.
	int r               = -1;
	const long double e = 1E-12L;
	const long double s = mc_raise2l(b) - 4.0L * a * c;

	if (mc_fabsl(s) < e) {
		*r1 = -0.5L * b / a;
		*i1 = 0.0L;
		*r2 = *r1;
		*i2 = *i1;
		 r  = 1;
	} else if (s > 0.0L) {
		const long double q = mc_sqrtl(s);
		*r1                 = (-b + q) * (1.0L / (a * 2.0L));
		*i1                 = 0.0L;
		*r2                 = (-b - q) * (1.0L / (a * 2.0L));
		*i2                 = 0.0L;
		 r                  = 2;
	} else if (s < 0.0L) {
		mc_zsqrtl(r1, i1, s, 0.0L);
		*r2 = *r1;
		*i2 = *i1;
		 r  = 4;

		mc_zaddl(r1, i1, -b, 0.0L, *r1, *i1);
		mc_zfdivl(r1, i1, *r1, *i1, a * 2.0);

		mc_zsubl(r2, i2, -b, 0.0L, *r2, *i2);
		mc_zfdivl(r2, i2, *r2, *i2, a * 2.0L);
	}
	return r;
}

#endif /* !MC_ZPOLYROOT2_H */

/* EOF */