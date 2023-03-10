//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_evalpoly.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fma.h>

#ifndef MC_EVALPOLY_H
#define MC_EVALPOLY_H

#pragma mark - mc_evalpoly2 -

MC_TARGET_PROC float mc_evalpoly2f(const float x
	, const float p1
	, const float p2
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	float s;
	if (f == 1) {
		s = p1;
		s = mc_fmaf(x, s, p2);
	} else {
		s = p2;
		s = mc_fmaf(x, s, p1);
	}
	return s;
}

MC_TARGET_PROC double mc_evalpoly2ff(const float x
	, const float p1
	, const float p2
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	const double xd  = mc_cast(double, x);
	const double p1d = mc_cast(double, p1);
	const double p2d = mc_cast(double, p2);

	double s;

	if (f == 1) {
		s = p1d;
		s = mc_fma(xd, s, p2d);
	} else {
		s = p2d;
		s = mc_fma(xd, s, p1d);
	}
	return s;
}

MC_TARGET_PROC double mc_evalpoly2(const double x
	, const double p1
	, const double p2
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	double s;
	if (f == 1) {
		s = p1;
		s = mc_fma(x, s, p2);
	} else {
		s = p2;
		s = mc_fma(x, s, p1);
	}
	return s;
}

MC_TARGET_PROC long double mc_evalpoly2l(const long double x
	, const long double p1
	, const long double p2
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	long double s;
	if (f == 1) {
		s = p1;
		s = mc_fmal(x, s, p2);
	} else {
		s = p2;
		s = mc_fmal(x, s, p1);
	}
	return s;
}

#pragma mark - mc_evalpoly3 -

MC_TARGET_PROC float mc_evalpoly3f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	float s;
	if (f == 1) {
		s = p1;
		s = mc_fmaf(x, s, p2);
		s = mc_fmaf(x, s, p3);
	} else {
		s = p3;
		s = mc_fmaf(x, s, p2);
		s = mc_fmaf(x, s, p1);
	}
	return s;
}

MC_TARGET_PROC double mc_evalpoly3ff(const float x
	, const float p1
	, const float p2
	, const float p3
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	const double xd  = mc_cast(double, x);
	const double p1d = mc_cast(double, p1);
	const double p2d = mc_cast(double, p2);
	const double p3d = mc_cast(double, p3);

	double s;

	if (f == 1) {
		s = p1d;
		s = mc_fma(xd, s, p2d);
		s = mc_fma(xd, s, p3d);
	} else {
		s = p3d;
		s = mc_fma(xd, s, p2d);
		s = mc_fma(xd, s, p1d);
	}
	return s;
}

MC_TARGET_PROC double mc_evalpoly3(const double x
	, const double p1
	, const double p2
	, const double p3
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	double s;
	if (f == 1) {
		s = p1;
		s = mc_fma(x, s, p2);
		s = mc_fma(x, s, p3);
	} else {
		s = p3;
		s = mc_fma(x, s, p2);
		s = mc_fma(x, s, p1);
	}
	return s;
}

MC_TARGET_PROC long double mc_evalpoly3l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	long double s;
	if (f == 1) {
		s = p1;
		s = mc_fmal(x, s, p2);
		s = mc_fmal(x, s, p3);
	} else {
		s = p3;
		s = mc_fmal(x, s, p2);
		s = mc_fmal(x, s, p1);
	}
	return s;
}

#pragma mark - mc_evalpoly4 -

MC_TARGET_PROC float mc_evalpoly4f(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	float s;
	if (f == 1) {
		s = p1;
		s = mc_fmaf(x, s, p2);
		s = mc_fmaf(x, s, p3);
		s = mc_fmaf(x, s, p4);
	} else {
		s = p4;
		s = mc_fmaf(x, s, p3);
		s = mc_fmaf(x, s, p2);
		s = mc_fmaf(x, s, p1);
	}
	return s;
}

MC_TARGET_PROC double mc_evalpoly4ff(const float x
	, const float p1
	, const float p2
	, const float p3
	, const float p4
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	const double xd  = mc_cast(double, x);
	const double p1d = mc_cast(double, p1);
	const double p2d = mc_cast(double, p2);
	const double p3d = mc_cast(double, p3);
	const double p4d = mc_cast(double, p4);

	double s;

	if (f == 1) {
		s = p1d;
		s = mc_fma(xd, s, p2d);
		s = mc_fma(xd, s, p3d);
		s = mc_fma(xd, s, p4d);
	} else {
		s = p4d;
		s = mc_fma(x, s, p3d);
		s = mc_fma(x, s, p2d);
		s = mc_fma(x, s, p1d);
	}
	return s;
}

MC_TARGET_PROC double mc_evalpoly4(const double x
	, const double p1
	, const double p2
	, const double p3
	, const double p4
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	double s;
	if (f == 1) {
		s = p1;
		s = mc_fma(x, s, p2);
		s = mc_fma(x, s, p3);
		s = mc_fma(x, s, p4);
	} else {
		s = p4;
		s = mc_fma(x, s, p3);
		s = mc_fma(x, s, p2);
		s = mc_fma(x, s, p1);
	}
	return s;
}

MC_TARGET_PROC long double mc_evalpoly4l(const long double x
	, const long double p1
	, const long double p2
	, const long double p3
	, const long double p4
	, const int f
) {
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	long double s;
	if (f == 1) {
		s = p1;
		s = mc_fmal(x, s, p2);
		s = mc_fmal(x, s, p3);
		s = mc_fmal(x, s, p4);
	} else {
		s = p4;
		s = mc_fmal(x, s, p3);
		s = mc_fmal(x, s, p2);
		s = mc_fmal(x, s, p1);
	}
	return s;
}

#pragma mark - mc_evalpoly -

MC_TARGET_PROC float mc_evalpolyf(const float x, const float * p, const int n, const int f)
{
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	float s = 0.0f;
	int i;
	if (mc_nonnullptr(p) && n > 0) {
		if (f == 1) {
			s = p[0];
			for(i = 1; i < n; i++) {
				s = mc_fmaf(x, s, p[i]);
			}
		} else {
			s = p[n - 1];
			for(i = n - 2; i >= 0; i--) {
				s = mc_fmaf(x, s, p[i]);
			}
		}
	}
	return s;
}

MC_TARGET_PROC double mc_evalpolyff(const float x, const float * p, const int n, const int f)
{
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	double s = 0.0;
	int i;
	if (mc_nonnullptr(p) && n > 0) {
		if (f == 1) {
			s = mc_cast(double, p[0]);
			for(i = 1; i < n; i++) {
				s = mc_fma(mc_cast(const double, x), s, mc_cast(const double, p[i]));
			}
		} else {
			s = mc_cast(double, p[n - 1]);
			for(i = n - 2; i >= 0; i--) {
				s = mc_fma(mc_cast(const double, x), s, mc_cast(const double, p[i]));
			}
		}
	}
	return s;
}

MC_TARGET_PROC double mc_evalpoly(const double x, const double * p, const int n, const int f)
{
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	double s = 0.0;
	int i;
	if (mc_nonnullptr(p) && n > 0) {
		if (f == 1) {
			s = p[0];
			for(i = 1; i < n; i++) {
				s = mc_fma(x, s, p[i]);
			}
		} else {
			s = p[n - 1];
			for(i = n - 2; i >= 0; i--) {
				s = mc_fma(x, s, p[i]);
			}
		}
	}
	return s;
}

MC_TARGET_PROC long double mc_evalpolyl(const long double x, const long double * p, const int n, const int f)
{
//!# Horner's method with `fma` computation.
//!# f=0: evaluating ascending power order.
//!# f=1: evaluating descending power order.
	long double s = 0.0L;
	int i;
	if (mc_nonnullptr(p) && n > 0) {
		if (f == 1) {
			s = p[0];
			for(i = 1; i < n; i++) {
				s = mc_fmal(x, s, p[i]);
			}
		} else {
			s = p[n - 1];
			for(i = n - 2; i >= 0; i--) {
				s = mc_fmal(x, s, p[i]);
			}
		}
	}
	return s;
}

#endif /* !MC_EVALPOLY_H */

/* EOF */