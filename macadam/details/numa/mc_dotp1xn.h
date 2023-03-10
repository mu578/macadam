//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_dotp1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fma.h>

#ifndef MC_DOTP1XN_H
#define MC_DOTP1XN_H

#pragma mark - mc_dotp1xn -

MC_TARGET_FUNC float mc_dotp1xnf(const int n, const float * x, const float * y, const int f)
{
#	if !MC_TARGET_HAVE_FMA
//!# TwoProduct split factor @see `mc_twoproduct`.
	const float cs = mc_cast_expr(const float, 4096 + 1);
#	endif

	int i   = 0;
	float s = 0.0f;

#	if !MC_TARGET_HAVE_FMA
	float w = 0.0f, h, q, r, x1, x2, y1, y2;
#	endif

	if (n > 0) {
		switch (f)
		{
			case 0:
				for (; i < n; i++) {
					s = s + (x[i] * y[i]);
				}
			break;
			case 1:
				for (; i < n; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fmaf(x[i], y[i], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(x[i],y[i],h,r).
					q  = x[i];
//!# split x[i] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = y[i];
//!# h=x[i]*y[i].
					h  = q * r;
//!# split y into y1,y2.
					q  = cs * r;
					y2 = q - r;
					y1 = q - y2;
					y2 = r - y1;
//!# r=x2*y2-(((h-x1*y1) - x2*y1) - x1*y2
					q  = x1 * y1;
					q  = h - q;
					y1 = y1 * x2;
					q  = q - y1;
					x1 = x1 * y2;
					q  = q - x1;
					x2 = x2 * y2;
					r  = x2 - q;
//!# (w,q)=TwoSum(w,h).
					x1 = w + h;
					x2 = x1 - w;
					y1 = x1 - x2;
					y2 = h - x2;
					q  = w - y1;
					q  = q + y2;
					w  = x1;
//!# s=s+(q+r).
					q = q + r;
					s = s + q;
#	endif
				}
			break;
		}
	}
#	if MC_TARGET_HAVE_FMA
	return s;
#	else
	return w + s;
#	endif
}

MC_TARGET_FUNC double mc_dotp1xnff(const int n, const float * x, const float * y, const int f)
{
#	if !MC_TARGET_HAVE_FMA
//!# TwoProduct split factor @see `mc_twoproduct`.
	const double cs = mc_cast_expr(const double, 134217728 + 1);
#	endif

	int i    = 0;
	double s = 0.0;

#	if !MC_TARGET_HAVE_FMA
	double w = 0.0, h, q, r, x1, x2, y1, y2;
#	endif

	if (n > 0) {
		switch (f)
		{
			case 0:
				for (; i < n; i++) {
					s = s + (mc_cast(double, x[i]) * mc_cast(double, y[i]));
				}
			break;
			case 1:
				for (; i < n; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fmafd(x[i], y[i], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(x[i],y[i],h,r).
					q  = mc_cast(double, x[i]);
//!# split x[i] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = mc_cast(double, y[i]);
//!# h=x[i]*y[i].
					h  = q * r;
//!# split y into y1,y2.
					q  = cs * r;
					y2 = q - r;
					y1 = q - y2;
					y2 = r - y1;
//!# r=x2*y2-(((h-x1*y1) - x2*y1) - x1*y2
					q  = x1 * y1;
					q  = h - q;
					y1 = y1 * x2;
					q  = q - y1;
					x1 = x1 * y2;
					q  = q - x1;
					x2 = x2 * y2;
					r  = x2 - q;
//!# (w,q)=TwoSum(w,h).
					x1 = w + h;
					x2 = x1 - w;
					y1 = x1 - x2;
					y2 = h - x2;
					q  = w - y1;
					q  = q + y2;
					w  = x1;
//!# s=s+(q+r).
					q = q + r;
					s = s + q;
#	endif
				}
			break;
		}
	}
#	if MC_TARGET_HAVE_FMA
	return s;
#	else
	return w + s;
#	endif
}

MC_TARGET_FUNC double mc_dotp1xn(const int n, const double * x, const double * y, const int f)
{
#	if !MC_TARGET_HAVE_FMA
//!# TwoProduct split factor @see `mc_twoproduct`.
	const double cs = mc_cast_expr(const double, 134217728 + 1);
#	endif

	int i    = 0;
	double s = 0.0;

#	if !MC_TARGET_HAVE_FMA
	double w = 0.0, h, q, r, x1, x2, y1, y2;
#	endif

	if (n > 0) {
		switch (f)
		{
			case 0:
				for (; i < n; i++) {
					s = s + (x[i] * y[i]);
				}
			break;
			case 1:
				for (; i < n; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fma(x[i], y[i], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(x[i],y[i],h,r).
					q  = x[i];
//!# split x[i] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = y[i];
//!# h=x[i]*y[i].
					h  = q * r;
//!# split y into y1,y2.
					q  = cs * r;
					y2 = q - r;
					y1 = q - y2;
					y2 = r - y1;
//!# r=x2*y2-(((h-x1*y1) - x2*y1) - x1*y2
					q  = x1 * y1;
					q  = h - q;
					y1 = y1 * x2;
					q  = q - y1;
					x1 = x1 * y2;
					q  = q - x1;
					x2 = x2 * y2;
					r  = x2 - q;
//!# (w,q)=TwoSum(w,h).
					x1 = w + h;
					x2 = x1 - w;
					y1 = x1 - x2;
					y2 = h - x2;
					q  = w - y1;
					q  = q + y2;
					w  = x1;
//!# s=s+(q+r).
					q = q + r;
					s = s + q;
#	endif
				}
			break;
		}
	}
#	if MC_TARGET_HAVE_FMA
	return s;
#	else
	return w + s;
#	endif
}

MC_TARGET_FUNC long double mc_dotp1xnl(const int n, const long double * x, const long double * y, const int f)
{
#	if !MC_TARGET_HAVE_FMA
//!# TwoProduct split factor @see `mc_twoproduct`.
#	if MC_TARGET_HAVE_LONG_DOUBLE && (MC_TARGET_LONG_DOUBLE_TYPE == MC_TARGET_LONG_DOUBLE_X87)
	const long double cs = mc_cast_expr(const long double, 4294967296 + 1);
#	elif MC_TARGET_HAVE_LONG_DOUBLE && (MC_TARGET_LONG_DOUBLE_TYPE == MC_TARGET_LONG_DOUBLE_ALIAS)
	const long double cs = mc_cast_expr(const long double, 134217728 + 1);
#	else
#	pragma message("Mantissa is too large. set @MC_TARGET_HAVE_FMA to 1.")
	const long double cs = MCK_NAN;
#	endif
#	endif

	int i         = 0;
	long double s = 0.0L;

#	if !MC_TARGET_HAVE_FMA
	long double  w = 0.0L, h, q, r, x1, x2, y1, y2;
#	endif

	if (n > 0) {
		switch (f)
		{
			case 0:
				for (; i < n; i++) {
					s = s + (x[i] * y[i]);
				}
			break;
			case 1:
				for (; i < n; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fmal(x[i], y[i], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(x[i],y[i],h,r).
					q  = x[i];
//!# split x[i] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = y[i];
//!# h=x[i]*y[i].
					h  = q * r;
//!# split y into y1,y2.
					q  = cs * r;
					y2 = q - r;
					y1 = q - y2;
					y2 = r - y1;
//!# r=x2*y2-(((h-x1*y1) - x2*y1) - x1*y2
					q  = x1 * y1;
					q  = h - q;
					y1 = y1 * x2;
					q  = q - y1;
					x1 = x1 * y2;
					q  = q - x1;
					x2 = x2 * y2;
					r  = x2 - q;
//!# (w,q)=TwoSum(w,h).
					x1 = w + h;
					x2 = x1 - w;
					y1 = x1 - x2;
					y2 = h - x2;
					q  = w - y1;
					q  = q + y2;
					w  = x1;
//!# s=s+(q+r).
					q = q + r;
					s = s + q;
#	endif
				}
			break;
		}
	}
#	if MC_TARGET_HAVE_FMA
	return s;
#	else
	return w + s;
#	endif
}

#endif /* !MC_DOTP1XN_H */

/* EOF */