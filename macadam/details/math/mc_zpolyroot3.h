//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zpolyroot3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_acos.h>
#include <macadam/details/math/mc_cbrt.h>
#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_raise3.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zsqrt.h>

#ifndef MC_ZPOLYROOT3_H
#define MC_ZPOLYROOT3_H

#pragma mark - mc_zpolyroot3 -

MC_TARGET_PROC int mc_zpolyroot3f(const float a, const float b, const float c, const float d
	, float * r1
	, float * i1
	, float * r2
	, float * i2
	, float * r3
	, float * i3
) {
//!# Three roots of the cubic equation: ax^3+bx^2+cx+cd=0.
//!# Returns 2 if roots are real and at least two are equal.
//!# Returns 3 if roots are real and distinct.
//!# Returns 5 if first root is real, the two others are conjugate.
	int r         = -1;
	float ca      = a;
	float cb      = b;
	float cc      = c;
	float cd      = d;
	const float e = 1E-05f;
	if (!(ca == 0.0f || cd == 0.0f)) {
		if (ca != 1.0f) {
			ca = 1.0f / ca;
			cb = cb * ca;
			cc = cc * ca;
			cd = cd * ca;
		}
		const float q = (3.0f * cc - mc_raise2f(cb)) / 9.0f;
		const float t = (9.0f * cb * cc - 27.0f * cd - 2.0f * mc_raise2f(cb) * cb) / 54.0f;
		const float s = mc_raise3f(q) + mc_raise2f(t);

		if (mc_fabsf(s) < e) {
			const float m = t > 0.0f ? mc_cbrtf(t) : -mc_cbrtf(-t);
			*r1           = 2.0f * m - cb / 3.0f;
			*i1           = 0.0f;
			*r2           = -m - cb / 3.0f;
			*i2           = 0.0f;
			*r3           = *r2;
			*i3           = *i2;
			 r            = 2;
		} else if (s < 0) {
			*r1 = -q;
			if (*r1 < 0.0f) {
				mc_zsqrtf(r1, i1, *r1, 0.0f);
			} else {
				*i1 = mc_sqrtf(*r1);
			}
			const float m = *i1;
			const float h = 2.0f * m;
			const float g = cb / 3.0f;
			const float w = mc_acosf(t / (m * mc_fabsf(q)));
			*r1           = h * mc_cosf(w / 3.0f) - g;
			*i1           = 0.0f;
			*r2           = h * mc_cosf((w + MCK_KF(MCK_2PI)) / 3.0f) - g;
			*i2           = *i1;
			*r3           = h * mc_cosf((w + MCK_KF(MCK_4PI)) / 3.0f) - g;
			*i3           = *i2;
			 r            = 3;
		} else if (s > 0) {
			float u, v, w, h, g;
			const float m = mc_sqrtf(s);
			 u            = t + m;
			 v            = t - m;
			 u            = (u > 0.0f) ? mc_cbrtf(u) : -mc_cbrtf(-u);
			 v            = (v > 0.0f) ? mc_cbrtf(v) : -mc_cbrtf(-v);
			 w            = u + v;
			 h            = MCK_KF(MCK_SQRT3) / 2.0f * (u - v);
			 g            = cb / 3.0f;
			*r1           = w - g;
			*i1           = 0.0f;
			*r2           = -0.5f * w - g;
			*i2           = h;
			*r3           = *r2;
			*i3           = -(*i2);
			 r            = 5;
		}
	}
	return r;
}

MC_TARGET_PROC int mc_zpolyroot3(const double a, const double b, const double c, const double d
	, double * r1
	, double * i1
	, double * r2
	, double * i2
	, double * r3
	, double * i3
) {
//!# Three roots of the cubic equation: ax^3+bx^2+cx+cd=0.
//!# Returns 2 if roots are real and at least two are equal.
//!# Returns 3 if roots are real and distinct.
//!# Returns 5 if first root is real, the two others are conjugate.
	int r          = -1;
	double ca      = a;
	double cb      = b;
	double cc      = c;
	double cd      = d;
	const double e = 1E-09;
	if (!(ca == 0.0 || cd == 0.0)) {
		if (ca != 1.0) {
			ca = 1.0 / ca;
			cb = cb * ca;
			cc = cc * ca;
			cd = cd * ca;
		}
		const double q = (3.0 * cc - mc_raise2(cb)) / 9.0;
		const double t = (9.0 * cb * cc - 27.0 * cd - 2.0 * mc_raise2(cb) * cb) / 54.0;
		const double s = mc_raise3(q) + mc_raise2(t);

		if (mc_fabs(s) < e) {
			const double m = t > 0.0 ? mc_cbrt(t) : -mc_cbrt(-t);
			*r1            = 2.0 * m - cb / 3.0;
			*i1            = 0.0;
			*r2            = -m - cb / 3.0;
			*i2            = 0.0;
			*r3            = *r2;
			*i3            = *i2;
			 r             = 2;
		} else if (s < 0) {
			*r1 = -q;
			if (*r1 < 0.0) {
				mc_zsqrt(r1, i1, *r1, 0.0);
			} else {
				*i1 = mc_sqrt(*r1);
			}
			const double m = *i1;
			const double h = 2.0 * m;
			const double g = cb / 3.0;
			const double w = mc_acos(t / (m * mc_fabs(q)));
			*r1            = h * mc_cos(w / 3.0) - g;
			*i1            = 0.0;
			*r2            = h * mc_cos((w + MCK_K(MCK_2PI)) / 3.0) - g;
			*i2            = *i1;
			*r3            = h * mc_cos((w + MCK_K(MCK_4PI)) / 3.0) - g;
			*i3            = *i2;
			 r             = 3;
		} else if (s > 0) {
			double u, v, w, h, g;
			const double m = mc_sqrt(s);
			 u             = t + m;
			 v             = t - m;
			 u             = (u > 0.0) ? mc_cbrt(u) : -mc_cbrt(-u);
			 v             = (v > 0.0) ? mc_cbrt(v) : -mc_cbrt(-v);
			 w             = u + v;
			 h             = MCK_K(MCK_SQRT3) / 2.0 * (u - v);
			 g             = cb / 3.0;
			*r1            = w - g;
			*i1            = 0.0;
			*r2            = -0.5 * w - g;
			*i2            = h;
			*r3            = *r2;
			*i3            = -(*i2);
			 r             = 5;
		}
	}
	return r;
}

MC_TARGET_PROC int mc_zpolyroot3l(const long double a, const long double b, const long double c, const long double d
	, long double * r1
	, long double * i1
	, long double * r2
	, long double * i2
	, long double * r3
	, long double * i3
) {
//!# Three roots of the cubic equation: ax^3+bx^2+cx+cd=0.
//!# Returns 2 if roots are real and at least two are equal.
//!# Returns 3 if roots are real and distinct.
//!# Returns 5 if first root is real, the two others are conjugate.
	int r               = -1;
	long double ca      = a;
	long double cb      = b;
	long double cc      = c;
	long double cd      = d;
	const long double e = 1E-12L;
	if (!(ca == 0.0L || cd == 0.0L)) {
		if (ca != 1.0L) {
			ca = 1.0L / ca;
			cb = cb * ca;
			cc = cc * ca;
			cd = cd * ca;
		}
		const long double q = (3.0L * cc - mc_raise2l(cb)) / 9.0L;
		const long double t = (9.0L * cb * cc - 27.0L * cd - 2.0L * mc_raise2l(cb) * cb) / 54.0L;
		const long double s = mc_raise3l(q) + mc_raise2l(t);

		if (mc_fabsl(s) < e) {
			const long double m = t > 0.0L ? mc_cbrtl(t) : -mc_cbrtl(-t);
			*r1                 = 2.0L * m - cb / 3.0L;
			*i1                 = 0.0L;
			*r2                 = -m - cb / 3.0L;
			*i2                 = 0.0L;
			*r3                 = *r2;
			*i3                 = *i2;
			 r                  = 2;
		} else if (s < 0) {
			*r1 = -q;
			if (*r1 < 0.0L) {
				mc_zsqrtl(r1, i1, *r1, 0.0L);
			} else {
				*i1 = mc_sqrtl(*r1);
			}
			const long double m = *i1;
			const long double h = 2.0L * m;
			const long double g = cb / 3.0L;
			const long double w = mc_acosl(t / (m * mc_fabsl(q)));
			*r1                 = h * mc_cosl(w / 3.0L) - g;
			*i1                 = 0.0L;
			*r2                 = h * mc_cosl((w + MCK_KL(MCK_2PI)) / 3.0L) - g;
			*i2                 = *i1;
			*r3                 = h * mc_cosl((w + MCK_KL(MCK_4PI)) / 3.0L) - g;
			*i3                 = *i2;
			 r                  = 3;
		} else if (s > 0) {
			long double u, v, w, h, g;
			const long double m = mc_sqrtl(s);
			 u                  = t + m;
			 v                  = t - m;
			 u                  = (u > 0.0L) ? mc_cbrtl(u) : -mc_cbrtl(-u);
			 v                  = (v > 0.0L) ? mc_cbrtl(v) : -mc_cbrtl(-v);
			 w                  = u + v;
			 h                  = MCK_KL(MCK_SQRT3) / 2.0L * (u - v);
			 g                  = cb / 3.0L;
			*r1                 = w - g;
			*i1                 = 0.0L;
			*r2                 = -0.5L * w - g;
			*i2                 = h;
			*r3                 = *r2;
			*i3                 = -(*i2);
			 r                  = 5;
		}
	}
	return r;
}

#endif /* !MC_ZPOLYROOT3_H */

/* EOF */