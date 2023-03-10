//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_summx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/numa/mc_2summx1.h>

#ifndef MC_SUMMX1_H
#define MC_SUMMX1_H

#pragma mark - mc_summx1 -

MC_TARGET_FUNC float mc_summx1f(const int m, const int n, const int j, const float * a, const int f)
{
	switch (f)
	{
		case 0:
		{
			int i   = 0;
			float s = 0.0f;
			for (; i < m; i++) {
				s = s + a[(n * i) + j];
			}
			return s;
		}
		case 1:
		{ //!# Kahan summation.
			int i   = 0;
			float s = 0.0f, c = 0.0f;
			for (; i < m; i++) {
				const float u = a[(n * i) + j] - c;
				const float w = s + u;
				c             = (w - s) - u;
				s             = w;
			}
			return s;
		}
		case 2:
		{ //!# Neumaier summation.
			int i   = 0;
			float s = 0.0f, c = 0.0f, u, w;
			for (; i < m; i++) {
				w = a[(n * i) + j];
				u = s + w;
				c = c + (mc_fabsf(s) >= mc_fabsf(w) ? (s - u) + w :  (w - u) + s);
				s = u;
			}
			return s + c;
		}
		case 3:
		{ //!# Klein summation.
			int i   = 0;
			float s = 0.0f, c = 0.0f, cs = 0.0f, ccs = 0.0f, cc, t, w;
			for (; i < m; i++) {
				w = a[(n * i) + j];
				t = s + w;
				if (mc_fabsf(s) >= mc_fabsf(w)) {
					c = (s - t) + w;
				} else {
					c = (w - t) + s;
				}
				s = t;
				t = cs + c;
				if (mc_fabsf(cs) >= mc_fabsf(c)) {
					cc = (cs - t) + c;
				} else {
					cc = (c - t) + cs;
				}
				cs  = t;
				ccs = ccs + cc;
			}
			return s + cs + ccs;
		}
		case 4:
		{//!# Absolute summation.
			int i   = 0;
			float s = 0.0f;
			for (; i < m; i++) {
				s = s + mc_fabsf(a[(n * i) + j]);
			}
			return s;
		}
		case 5:
		{//!# 2sum
			return mc_2summx1f(m, n, j, a);
		}
	}
	return 0.0f;
}

MC_TARGET_FUNC double mc_summx1ff(const int m, const int n, const int j, const float * a, const int f)
{
	switch (f)
	{
		case 0:
		{
			int i    = 0;
			double s = 0.0;
			for (; i < m; i++) {
				s = s + mc_cast(float, a[(n * i) + j]);
			}
			return s;
		}
		case 1:
		{ //!# Kahan summation.
			int i    = 0;
			double s = 0.0, c = 0.0;
			for (; i < m; i++) {
				const double u = mc_cast(double, a[(n * i) + j]) - c;
				const double w = s + u;
				c              = (w - s) - u;
				s              = w;
			}
			return s;
		}
		case 2:
		{ //!# Neumaier summation.
			int i    = 0;
			double s = 0.0, c = 0.0, u, w;
			for (; i < m; i++) {
				w = mc_cast(double, a[(n * i) + j]);
				u = s + w;
				c = c + (mc_fabs(s) >= mc_fabs(w) ? (s - u) + w :  (w - u) + s);
				s = u;
			}
			return s + c;
		}
		case 3:
		{ //!# Klein summation.
			int i    = 0;
			double s = 0.0, c = 0.0, cs = 0.0, ccs = 0.0, cc, t, w;
			for (; i < m; i++) {
				w = mc_cast(double, a[(n * i) + j]);
				t = s + w;
				if (mc_fabs(s) >= mc_fabs(w)) {
					c = (s - t) + w;
				} else {
					c = (w - t) + s;
				}
				s = t;
				t = cs + c;
				if (mc_fabs(cs) >= mc_fabs(c)) {
					cc = (cs - t) + c;
				} else {
					cc = (c - t) + cs;
				}
				cs  = t;
				ccs = ccs + cc;
			}
			return s + cs + ccs;
		}
		case 4:
		{//!# Absolute summation.
			int i    = 0;
			double s = 0.0;
			for (; i < m; i++) {
				s = s + mc_fabs(mc_cast(float, a[(n * i) + j]));
			}
			return s;
		}
		case 5:
		{//!# 2sum
			return mc_2summx1ff(m, n, j, a);
		}
	}
	return 0.0;
}

MC_TARGET_FUNC double mc_summx1(const int m, const int n, const int j, const double * a, const int f)
{
	switch (f)
	{
		case 0:
		{
			int i    = 0;
			double s = 0.0;
			for (; i < m; i++) {
				s = s + a[(n * i) + j];
			}
			return s;
		}
		case 1:
		{ //!# Kahan summation.
			int i    = 0;
			double s = 0.0, c = 0.0;
			for (; i < m; i++) {
				const double u = a[(n * i) + j] - c;
				const double w = s + u;
				c              = (w - s) - u;
				s              = w;
			}
			return s;
		}
		case 2:
		{ //!# Neumaier summation.
			int i    = 0;
			double s = 0.0, c = 0.0, u, w;
			for (; i < m; i++) {
				w = a[(n * i) + j];
				u = s + w;
				c = c + (mc_fabs(s) >= mc_fabs(w) ? (s - u) + w :  (w - u) + s);
				s = u;
			}
			return s + c;
		}
		case 3:
		{ //!# Klein summation.
			int i    = 0;
			double s = 0.0, c = 0.0, cs = 0.0, ccs = 0.0, cc, t, w;
			for (; i < m; i++) {
				w = a[(n * i) + j];
				t = s + w;
				if (mc_fabs(s) >= mc_fabs(w)) {
					c = (s - t) + w;
				} else {
					c = (w - t) + s;
				}
				s = t;
				t = cs + c;
				if (mc_fabs(cs) >= mc_fabs(c)) {
					cc = (cs - t) + c;
				} else {
					cc = (c - t) + cs;
				}
				cs  = t;
				ccs = ccs + cc;
			}
			return s + cs + ccs;
		}
		case 4:
		{//!# Absolute summation.
			int i    = 0;
			double s = 0.0;
			for (; i < m; i++) {
				s = s + mc_fabs(a[(n * i) + j]);
			}
			return s;
		}
		case 5:
		{//!# 2sum
			return mc_2summx1(m, n, j, a);
		}
	}
	return 0.0;
}

MC_TARGET_FUNC long double mc_summx1l(const int m, const int n, const int j, const long double * a, const int f)
{
	switch (f)
	{
		case 0:
		{
			int i         = 0;
			long double s = 0.0L;
			for (; i < m; i++) {
				s = s + a[(n * i) + j];
			}
			return s;
		}
		case 1:
		{ //!# Kahan summation.
			int i         = 0;
			long double s = 0.0L, c = 0.0L;
			for (; i < m; i++) {
				const long double u = a[(n * i) + j] - c;
				const long double w = s + u;
				c                   = (w - s) - u;
				s                   = w;
			}
			return s;
		}
		case 2:
		{ //!# Neumaier summation.
			int i         = 0;
			long double s = 0.0L, c = 0.0L, u, w;
			for (; i < m; i++) {
				w = a[(n * i) + j];
				u = s + w;
				c = c + (mc_fabsl(s) >= mc_fabsl(w) ? (s - u) + w :  (w - u) + s);
				s = u;
			}
			return s + c;
		}
		case 3:
		{ //!# Klein summation.
			int i         = 0;
			long double s = 0.0L, c = 0.0L, cs = 0.0L, ccs = 0.0L, cc, t, w;
			for (; i < m; i++) {
				w = a[(n * i) + j];
				t = s + w;
				if (mc_fabsl(s) >= mc_fabsl(w)) {
					c = (s - t) + w;
				} else {
					c = (w - t) + s;
				}
				s = t;
				t = cs + c;
				if (mc_fabsl(cs) >= mc_fabsl(c)) {
					cc = (cs - t) + c;
				} else {
					cc = (c - t) + cs;
				}
				cs  = t;
				ccs = ccs + cc;
			}
			return s + cs + ccs;
		}
		case 4:
		{//!# Absolute summation.
			int i         = 0;
			long double s = 0.0L;
			for (; i < m; i++) {
				s = s + mc_fabsl(a[(n * i) + j]);
			}
			return s;
		}
		case 5:
		{//!# 2sum
			return mc_2summx1l(m, n, j, a);
		}
	}
	return 0.0L;
}

#endif /* !MC_SUMMX1_H */

/* EOF */