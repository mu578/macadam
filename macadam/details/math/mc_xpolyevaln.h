//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_xpolyevaln.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>

#ifndef MC_XPOLYEVALN_H
#define MC_XPOLYEVALN_H

#pragma mark - mc_xpolyevalne -

MC_TARGET_PROC float mc_xpolyevalnef(const float x, const float * p, const unsigned int n, float * err)
{
//!# Evaluating the polynomial p (in ascending powers order) at
//!# the specified value of x + computes absolute error estimate.
	int i, m;
	float z, s = 0.0f, e;
	if (n > 0 && n < 256) {
		m = mc_cast(int, n);
		s = p[m - 1];
		if (mc_nonnullptr(err)) {
			z = mc_fabsf(x);
			e = 0.5f * mc_fabsf(s);
			for (i = m - 2; i >= 0; i--) {
				s = s * x + p[i];
				e = e * z + mc_fabsf(s);
			}
			 e   = MCLIMITS_EPSILONF * mc_fabsf(2.0f * e - mc_fabsf(s));
			*err = e;
		} else {
			for (i = m - 2; i >= 0; i--) {
				s = s * x + p[i];
			}
		}
	}
	return s;
}

MC_TARGET_PROC double mc_xpolyevalne(const double x, const double * p, const unsigned int n, double * err)
{
//!# Evaluating the polynomial p (in ascending powers order) at
//!# the specified value of x + computes absolute error estimate.
	int i, m;
	double z, s = 0.0, e;
	if (n > 0 && n < 256) {
		m = mc_cast(int, n);
		s = p[m - 1];
		if (mc_nonnullptr(err)) {
			z = mc_fabs(x);
			e = 0.5 * mc_fabs(s);
			for (i = m - 2; i >= 0; i--) {
				s = s * x + p[i];
				e = e * z + mc_fabs(s);
			}
			 e   = MCLIMITS_EPSILON * mc_fabs(2.0 * e - mc_fabs(s));
			*err = e;
		} else {
			for (i = m - 2; i >= 0; i--) {
				s = s * x + p[i];
			}
		}
	}
	return s;
}

MC_TARGET_PROC long double mc_xpolyevalnel(const long double x, const long double * p, const unsigned int n, long double * err)
{
//!# Evaluating the polynomial p (in ascending powers order) at
//!# the specified value of x + computes absolute error estimate.
	int i, m;
	long double z, s = 0.0L, e;
	if (n > 0 && n < 256) {
		m = mc_cast(int, n);
		s = p[m - 1];
		if (mc_nonnullptr(err)) {
			z = mc_fabsl(x);
			e = 0.5L * mc_fabsl(s);
			for (i = m - 2; i >= 0; i--) {
				s = s * x + p[i];
				e = e * z + mc_fabsl(s);
			}
			 e   = MCLIMITS_EPSILONL * mc_fabsl(2.0L * e - mc_fabsl(s));
			*err = e;
		} else {
			for (i = m - 2; i >= 0; i--) {
				s = s * x + p[i];
			}
		}
	}
	return s;
}

#pragma mark - mc_xpolyeval2 -

MC_TARGET_PROC float mc_xpolyeval2f(const float x
	, const float p1
	, const float p2
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval2(const double x
	, const double p1
	, const double p2
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval2l(const long double x
	, const long double p1
	, const long double p2
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval3 -

MC_TARGET_PROC float mc_xpolyeval3f(const float x
	, const float p1
	, const float p2
	, const float p3
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval3(const double x
	, const double p1
	, const double p2
	, const double p3
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval3l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval4 -

MC_TARGET_PROC float mc_xpolyeval4f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval4(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval4l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval5 -

MC_TARGET_PROC float mc_xpolyeval5f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval5(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval5l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval6 -

MC_TARGET_PROC float mc_xpolyeval6f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval6(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval6l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval7 -

MC_TARGET_PROC float mc_xpolyeval7f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval7(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval7l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval8 -

MC_TARGET_PROC float mc_xpolyeval8f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval8(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval8l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval9 -

MC_TARGET_PROC float mc_xpolyeval9f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval9(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval9l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval10 -

MC_TARGET_PROC float mc_xpolyeval10f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval10(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval10l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval11 -

MC_TARGET_PROC float mc_xpolyeval11f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval11(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval11l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval12 -

MC_TARGET_PROC float mc_xpolyeval12f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval12(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval12l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval13 -

MC_TARGET_PROC float mc_xpolyeval13f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval13(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval13l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval14 -

MC_TARGET_PROC float mc_xpolyeval14f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
	, const float p14
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval14(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
	, const double p14
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval14l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
	, const long double p14
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval15 -

MC_TARGET_PROC float mc_xpolyeval15f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
	, const float p14
	, const float p15
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval15(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
	, const double p14
	, const double p15
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval15l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
	, const long double p14
	, const long double p15
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval16 -

MC_TARGET_PROC float mc_xpolyeval16f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
	, const float p14
	, const float p15
	, const float p16
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval16(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
	, const double p14
	, const double p15
	, const double p16
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval16l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
	, const long double p14
	, const long double p15
	, const long double p16
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval17 -

MC_TARGET_PROC float mc_xpolyeval17f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
	, const float p14
	, const float p15
	, const float p16
	, const float p17
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval17(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
	, const double p14
	, const double p15
	, const double p16
	, const double p17
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval17l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
	, const long double p14
	, const long double p15
	, const long double p16
	, const long double p17
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval18 -

MC_TARGET_PROC float mc_xpolyeval18f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
	, const float p14
	, const float p15
	, const float p16
	, const float p17
	, const float p18
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval18(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
	, const double p14
	, const double p15
	, const double p16
	, const double p17
	, const double p18
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval18l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
	, const long double p14
	, const long double p15
	, const long double p16
	, const long double p17
	, const long double p18
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval19 -

MC_TARGET_PROC float mc_xpolyeval19f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
	, const float p14
	, const float p15
	, const float p16
	, const float p17
	, const float p18
	, const float p19
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p19;
	s = s * x + p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval19(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
	, const double p14
	, const double p15
	, const double p16
	, const double p17
	, const double p18
	, const double p19
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p19;
	s = s * x + p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval19l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
	, const long double p14
	, const long double p15
	, const long double p16
	, const long double p17
	, const long double p18
	, const long double p19
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p19;
	s = s * x + p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyeval20 -

MC_TARGET_PROC float mc_xpolyeval20f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const float p5
	, const float p6
	, const float p7
	, const float p8
	, const float p9
	, const float p10
	, const float p11
	, const float p12
	, const float p13
	, const float p14
	, const float p15
	, const float p16
	, const float p17
	, const float p18
	, const float p19
	, const float p20
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	s = p20;
	s = s * x + p19;
	s = s * x + p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC double mc_xpolyeval20(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const double p5
	, const double p6
	, const double p7
	, const double p8
	, const double p9
	, const double p10
	, const double p11
	, const double p12
	, const double p13
	, const double p14
	, const double p15
	, const double p16
	, const double p17
	, const double p18
	, const double p19
	, const double p20
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	s = p20;
	s = s * x + p19;
	s = s * x + p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

MC_TARGET_PROC long double mc_xpolyeval20l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const long double p5
	, const long double p6
	, const long double p7
	, const long double p8
	, const long double p9
	, const long double p10
	, const long double p11
	, const long double p12
	, const long double p13
	, const long double p14
	, const long double p15
	, const long double p16
	, const long double p17
	, const long double p18
	, const long double p19
	, const long double p20
) {
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	s = p20;
	s = s * x + p19;
	s = s * x + p18;
	s = s * x + p17;
	s = s * x + p16;
	s = s * x + p15;
	s = s * x + p14;
	s = s * x + p13;
	s = s * x + p12;
	s = s * x + p11;
	s = s * x + p10;
	s = s * x + p9;
	s = s * x + p8;
	s = s * x + p7;
	s = s * x + p6;
	s = s * x + p5;
	s = s * x + p4;
	s = s * x + p3;
	s = s * x + p2;
	s = s * x + p1;
	return s;
}

#pragma mark - mc_xpolyevaln -

MC_TARGET_PROC float mc_xpolyevalnf(const float x, const float * p, const unsigned int n)
{
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	float s;
	switch (n)
	{
	case 2:
		s = mc_xpolyeval2f(x, p[0], p[1]);
	break;
	case 3:
		s = mc_xpolyeval3f(x, p[0], p[1], p[2]);
	break;
	case 4:
		s = mc_xpolyeval4f(x, p[0], p[1], p[2], p[3]);
	break;
	case 5:
		s = mc_xpolyeval5f(x, p[0], p[1], p[2], p[3], p[4]);
	break;
	case 6:
		s = mc_xpolyeval6f(x, p[0], p[1], p[2], p[3], p[4], p[5]);
	break;
	case 7:
		s = mc_xpolyeval7f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
	break;
	case 8:
		s = mc_xpolyeval8f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
	break;
	case 9:
		s = mc_xpolyeval9f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
	break;
	case 10:
		s = mc_xpolyeval10f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]);
	break;
	case 11:
		s = mc_xpolyeval11f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]);
	break;
	case 12:
		s = mc_xpolyeval12f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11]);
	break;
	case 13:
		s = mc_xpolyeval13f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12]);
	break;
	case 14:
		s = mc_xpolyeval14f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13]);
	break;
	case 15:
		s = mc_xpolyeval15f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14]);
	break;
	case 16:
		s = mc_xpolyeval16f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]);
	break;
	case 17:
		s = mc_xpolyeval17f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16]);
	break;
	case 18:
		s = mc_xpolyeval18f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17]);
	break;
	case 19:
		s = mc_xpolyeval19f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18]);
	break;
	case 20:
		s = mc_xpolyeval20f(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18], p[19]);
	break;
	default:
		s = mc_xpolyevalnef(x, p, n, MC_NULLPTR);
	}
	return s;
}

MC_TARGET_PROC double mc_xpolyevaln(const double x, const double * p, const unsigned int n)
{
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	double s;
	switch (n)
	{
	case 2:
		s = mc_xpolyeval2(x, p[0], p[1]);
	break;
	case 3:
		s = mc_xpolyeval3(x, p[0], p[1], p[2]);
	break;
	case 4:
		s = mc_xpolyeval4(x, p[0], p[1], p[2], p[3]);
	break;
	case 5:
		s = mc_xpolyeval5(x, p[0], p[1], p[2], p[3], p[4]);
	break;
	case 6:
		s = mc_xpolyeval6(x, p[0], p[1], p[2], p[3], p[4], p[5]);
	break;
	case 7:
		s = mc_xpolyeval7(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
	break;
	case 8:
		s = mc_xpolyeval8(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
	break;
	case 9:
		s = mc_xpolyeval9(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
	break;
	case 10:
		s = mc_xpolyeval10(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]);
	break;
	case 11:
		s = mc_xpolyeval11(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]);
	break;
	case 12:
		s = mc_xpolyeval12(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11]);
	break;
	case 13:
		s = mc_xpolyeval13(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12]);
	break;
	case 14:
		s = mc_xpolyeval14(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13]);
	break;
	case 15:
		s = mc_xpolyeval15(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14]);
	break;
	case 16:
		s = mc_xpolyeval16(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]);
	break;
	case 17:
		s = mc_xpolyeval17(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[17]);
	break;
	case 18:
		s = mc_xpolyeval18(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17]);
	break;
	case 19:
		s = mc_xpolyeval19(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18]);
	break;
	case 20:
		s = mc_xpolyeval20(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18], p[19]);
	break;
	default:
		s = mc_xpolyevalne(x, p, n, MC_NULLPTR);
	}
	return s;
}

MC_TARGET_PROC long double mc_xpolyevalnl(const long double x, const long double * p, const unsigned int n)
{
//!# Evaluating the polynomial p (in ascending powers order) at the specified value of x.
	long double s;
	switch (n)
	{
	case 2:
		s = mc_xpolyeval2l(x, p[0], p[1]);
	break;
	case 3:
		s = mc_xpolyeval3l(x, p[0], p[1], p[2]);
	break;
	case 4:
		s = mc_xpolyeval4l(x, p[0], p[1], p[2], p[3]);
	break;
	case 5:
		s = mc_xpolyeval5l(x, p[0], p[1], p[2], p[3], p[4]);
	break;
	case 6:
		s = mc_xpolyeval6l(x, p[0], p[1], p[2], p[3], p[4], p[5]);
	break;
	case 7:
		s = mc_xpolyeval7l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
	break;
	case 8:
		s = mc_xpolyeval8l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
	break;
	case 9:
		s = mc_xpolyeval9l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
	break;
	case 10:
		s = mc_xpolyeval10l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]);
	break;
	case 11:
		s = mc_xpolyeval11l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]);
	break;
	case 12:
		s = mc_xpolyeval12l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11]);
	break;
	case 13:
		s = mc_xpolyeval13l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12]);
	break;
	case 14:
		s = mc_xpolyeval14l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13]);
	break;
	case 15:
		s = mc_xpolyeval15l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14]);
	break;
	case 16:
		s = mc_xpolyeval16l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]);
	break;
	case 17:
		s = mc_xpolyeval17l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16]);
	break;
	case 18:
		s = mc_xpolyeval18l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17]);
	break;
	case 19:
		s = mc_xpolyeval19l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18]);
	break;
	case 20:
		s = mc_xpolyeval20l(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18], p[19]);
	break;
	default:
		s = mc_xpolyevalnel(x, p, n, MC_NULLPTR);
	}
	return s;
}

#endif /* !MC_XPOLYEVALN_H */

/* EOF */