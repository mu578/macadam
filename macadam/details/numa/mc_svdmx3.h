//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_svdmx3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_svdgr1mxn.h>
#include <macadam/details/numa/mc_zerosnxn.h>

#ifndef MC_SVDMX3_H
#define MC_SVDMX3_H

#pragma mark - mc_svdmx3 -

MC_TARGET_FUNC int mc_svdmx3f(const int m, const float * a, float * MC_TARGET_RESTRICT u, float s[9], float v[9])
{
	const float eps = (m / 10.0f) * mc_sqrtf(MCLIMITS_EPSILONF);
	const float tol = 10E-08f;
	float * w       = s + 3;
	if (m > 2 && 0 == mc_svdgr1mxnf(m, 3, a, w, 1, 1, 1, eps, tol, u, s, v)) {
		const float sv1 = s[0];
		const float sv2 = s[1];
		const float sv3 = s[2];
		mc_zerosnxnf(3, s);
		s[0]            = sv1;
		s[4]            = sv2;
		s[8]            = sv3;
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svdmx3ff(const int m, const float * a, double * MC_TARGET_RESTRICT u, double s[9], double v[9])
{
	const double eps = (m / 10.0) * mc_sqrt(MCLIMITS_EPSILON);
	const double tol = 10E-12;
	double * w       = s + 3;
	if (m > 2 && 0 == mc_svdgr1mxnff(m, 3, a, w, 1, 1, 1, eps, tol, u, s, v)) {
		const double sv1 = s[0];
		const double sv2 = s[1];
		const double sv3 = s[2];
		mc_zerosnxn(3, s);
		s[0]             = sv1;
		s[4]             = sv2;
		s[8]             = sv3;
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svdmx3(const int m, const double * a, double * MC_TARGET_RESTRICT u, double s[9], double v[9])
{
	const double eps = (m / 10.0) * mc_sqrt(MCLIMITS_EPSILON);
	const double tol = 10E-12;
	double * w       = s + 3;
	if (m > 2 && 0 == mc_svdgr1mxn(m, 3, a, w, 1, 1, 1, eps, tol, u, s, v)) {
		const double sv1 = s[0];
		const double sv2 = s[1];
		const double sv3 = s[2];
		mc_zerosnxn(3, s);
		s[0]             = sv1;
		s[4]             = sv2;
		s[8]             = sv3;
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svdmx3l(const int m, const long double * a, long double * MC_TARGET_RESTRICT u, long double s[9], long double v[9])
{
	const long double eps = (m / 10.0L) * mc_sqrtl(MCLIMITS_EPSILONL);
	const long double tol = 10E-15L;
	long double * w       = s + 3;
	if (m > 2 && 0 == mc_svdgr1mxnl(m, 3, a, w, 1, 1, 1, eps, tol, u, s, v)) {
		const long double sv1 = s[0];
		const long double sv2 = s[1];
		const long double sv3 = s[2];
		mc_zerosnxnl(3, s);
		s[0]                  = sv1;
		s[4]                  = sv2;
		s[8]                  = sv3;
		return 0;
	}
	return -1;
}

#endif /* !MC_SVDMX3_H */

/* EOF */