//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_chpoly2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_det2x2.h>
#include <macadam/details/numa/mc_trace2x2.h>

#ifndef MC_CHPOLY2X2_H
#define MC_CHPOLY2X2_H

#pragma mark - mc_chpoly2x2 -

MC_TARGET_FUNC void mc_chpoly2x2f(const float a[3], float * p0, float * p1, float * p2)
{
//!#
//!# cp=($p0*lambda^2) + ($p1*lambda) + ($p2)
//!# $p0=1, $p1=trace(a), $p2=det(a)
//!#
	*p0 = 1.0f;
	*p1 = mc_trace2x2f(a);
	*p2 = mc_det2x2f(a);
}

MC_TARGET_FUNC void mc_chpoly2x2ff(const float a[3], double * p0, double * p1, double * p2)
{
//!#
//!# cp=($p0*lambda^2) + ($p1*lambda) + ($p2)
//!# $p0=1, $p1=trace(a), $p2=det(a)
//!#
	*p0 = 1.0;
	*p1 = mc_trace2x2ff(a);
	*p2 = mc_det2x2ff(a);
}

MC_TARGET_FUNC void mc_chpoly2x2(const double a[3], double * p0, double * p1, double * p2)
{
//!#
//!# cp=($p0*lambda^2) + ($p1*lambda) + ($p2)
//!# $p0=1, $p1=trace(a), $p2=det(a)
//!#
	*p0 = 1.0;
	*p1 = mc_trace2x2(a);
	*p2 = mc_det2x2(a);
}

MC_TARGET_FUNC void mc_chpoly2x2l(const long double a[3], long double * p0, long double * p1, long double * p2)
{
//!#
//!# cp=($p0*lambda^2) + ($p1*lambda) + ($p2)
//!# $p0=1, $p1=trace(a), $p2=det(a)
//!#
	*p0 = 1.0L;
	*p1 = mc_trace2x2l(a);
	*p2 = mc_det2x2l(a);
}

#endif /* !MC_CHPOLY2X2_H */

/* EOF */