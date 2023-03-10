//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_chpoly3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/numa/mc_det3x3.h>
#include <macadam/details/numa/mc_trace3x3.h>

#ifndef MC_CHPOLY3X3_H
#define MC_CHPOLY3X3_H

#pragma mark - mc_chpoly3x3 -

MC_TARGET_FUNC void mc_chpoly3x3f(const float a[3], float * p0, float * p1, float * p2, float * p3)
{
//!#
//!# cp=($p0*lambda^3) + ($p1*lambda^2) + ($p2*lambda) + ($p3)
//!# $p0=-1, $p1=trace(a), $p2=$p1^2 - trace(a^2), $p3=det(a)
//!#
	*p0 = -1.0f;
	*p1 = mc_trace3x3f(a);
	*p2 = mc_raise2f(*p1) - (
		  (a[0] * a[0]) + (a[1] * a[3]) + (a[2] * a[6])
		+ (a[3] * a[1]) + (a[4] * a[4]) + (a[5] * a[7])
		+ (a[6] * a[2]) + (a[7] * a[5]) + (a[8] * a[8])
	);
	*p3 = mc_det3x3f(a);
}

MC_TARGET_FUNC void mc_chpoly3x3ff(const float a[3], double * p0, double * p1, double * p2, double * p3)
{
//!#
//!# cp=($p0*lambda^3) + ($p1*lambda^2) + ($p2*lambda) + ($p3)
//!# $p0=-1, $p1=trace(a), $p2=$p1^2 - trace(a^2), $p3=det(a)
//!#
	*p0 = -1.0;
	*p1 = mc_trace3x3ff(a);
	*p2 = mc_raise2(*p1) - (
		  (mc_cast(double, a[0]) * mc_cast(double, a[0])) + (mc_cast(double, a[1]) * mc_cast(double, a[3])) + (mc_cast(double, a[2]) * mc_cast(double, a[6]))
		+ (mc_cast(double, a[3]) * mc_cast(double, a[1])) + (mc_cast(double, a[4]) * mc_cast(double, a[4])) + (mc_cast(double, a[5]) * mc_cast(double, a[7]))
		+ (mc_cast(double, a[6]) * mc_cast(double, a[2])) + (mc_cast(double, a[7]) * mc_cast(double, a[5])) + (mc_cast(double, a[8]) * mc_cast(double, a[8]))
	);
	*p3 = mc_det3x3ff(a);
}

MC_TARGET_FUNC void mc_chpoly3x3(const double a[3], double * p0, double * p1, double * p2, double * p3)
{
//!#
//!# cp=($p0*lambda^3) + ($p1*lambda^2) + ($p2*lambda) + ($p3)
//!# $p0=-1, $p1=trace(a), $p2=$p1^2 - trace(a^2), $p3=det(a)
//!#
	*p0 = -1.0;
	*p1 = mc_trace3x3(a);
	*p2 = mc_raise2(*p1) - (
		  (a[0] * a[0]) + (a[1] * a[3]) + (a[2] * a[6])
		+ (a[3] * a[1]) + (a[4] * a[4]) + (a[5] * a[7])
		+ (a[6] * a[2]) + (a[7] * a[5]) + (a[8] * a[8])
	);
	*p3 = mc_det3x3(a);
}

MC_TARGET_FUNC void mc_chpoly3x3l(const long double a[3], long double * p0, long double * p1, long double * p2, long double * p3)
{
//!#
//!# cp=($p0*lambda^3) + ($p1*lambda^2) + ($p2*lambda) + ($p3)
//!# $p0=-1, $p1=trace(a), $p2=$p1^2 - trace(a^2), $p3=det(a)
//!#
	*p0 = -1.0L;
	*p1 = mc_trace3x3l(a);
	*p2 = mc_raise2l(*p1) - (
		  (a[0] * a[0]) + (a[1] * a[3]) + (a[2] * a[6])
		+ (a[3] * a[1]) + (a[4] * a[4]) + (a[5] * a[7])
		+ (a[6] * a[2]) + (a[7] * a[5]) + (a[8] * a[8])
	);
	*p3 = mc_det3x3l(a);
}

#endif /* !MC_CHPOLY3X3_H */

/* EOF */