//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_dotpmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_DOTPMX1_H
#define MC_DOTPMX1_H

#pragma mark - mc_dotpmx1 -

MC_TARGET_FUNC float mc_dotpmx1f(const int m, const int n, const int p, const int j, const int k, const float * a, const float * b, const int f)
{
//!# Requires a[m x n] and b[m x p].
//!# A and B may be the same.
#	if !MC_TARGET_HAVE_FMA
//!# TwoProduct split factor @see `mc_twoproduct`.
	const float cs = mc_cast_expr(float, 4096 + 1);
#	endif

	int i   = 0;
	float s = 0.0f;

#	if !MC_TARGET_HAVE_FMA
	float w = 0.0f, h, q, r, x1, x2, y1, y2;
#	endif

	if (m > 0)
	{
		switch (f) {
			case 0:
				for (; i < m; i++) {
					s = s + (a[(n * i) + j] * b[(p * i) + k]);
				}
			break;
			case 1:
				for (; i < m; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fmaf(a[(n * i) + j], b[(p * i) + k], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(a[(n * i) + j],b[(p * i) + k],h,r).
					q  = a[(n * i) + j];
//!# split a[(n * i) + j] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = b[(p * i) + k];
//!# h=a[(n * i) + j]*b[(p * i) + k].
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

MC_TARGET_FUNC double mc_dotpmx1ff(const int m, const int n, const int p, const int j, const int k, const float * a, const float * b, const int f)
{
//!# Requires a[m x n] and b[m x p].
//!# A and B may be the same.
#	if !MC_TARGET_HAVE_FMA
//!# TwoProduct split factor @see `mc_twoproduct`.
	const double cs = mc_cast_expr(const double, 134217728 + 1);
#	endif

	int i    = 0;
	double s = 0.0;

#	if !MC_TARGET_HAVE_FMA
	double w = 0.0, h, q, r, x1, x2, y1, y2;
#	endif

	if (m > 0) {
		switch (f)
		{
			case 0:
				for (; i < m; i++) {
					s = s + (mc_cast(double, a[(n * i) + j]) * mc_cast(double, b[(p * i) + k]));
				}
			break;
			case 1:
				for (; i < m; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fmafd(a[(n * i) + j], b[(p * i) + k], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(a[(n * i) + j],b[(p * i) + k],h,r).
					q  = mc_cast(double, a[(n * i) + j]);
//!# split a[(n * i) + j] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = mc_cast(double, b[(p * i) + k]);
//!# h=a[(n * i) + j]*b[(p * i) + k].
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

MC_TARGET_FUNC double mc_dotpmx1(const int m, const int n, const int p, const int j, const int k, const double * a, const double * b, const int f)
{
//!# Requires a[m x n] and b[m x p].
//!# A and B may be the same.
#	if !MC_TARGET_HAVE_FMA
//!# TwoProduct split factor @see `mc_twoproduct`.
	const double cs = mc_cast_expr(const double, 134217728 + 1);
#	endif

	int i    = 0;
	double s = 0.0;

#	if !MC_TARGET_HAVE_FMA
	double w = 0.0, h, q, r, x1, x2, y1, y2;
#	endif

	if (m > 0) {
		switch (f)
		{
			case 0:
				for (; i < m; i++) {
					s = s + (a[(n * i) + j] * b[(p * i) + k]);
				}
			break;
			case 1:
				for (; i < m; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fma(a[(n * i) + j], b[(p * i) + k], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(a[(n * i) + j],b[(p * i) + k],h,r).
					q  = a[(n * i) + j];
//!# split a[(n * i) + j] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = b[(p * i) + k];
//!# h=a[(n * i) + j]*b[(p * i) + k].
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

MC_TARGET_FUNC long double mc_dotpmx1l(const int m, const int n, const int p, const int j, const int k, const long double * a, const long double * b, const int f)
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
	long double w = 0.0L, h, q, r, x1, x2, y1, y2;
#	endif

	if (m > 0) {
		switch (f)
		{
			case 0:
				for (; i < m; i++) {
					s = s + (a[(n * i) + j] * b[(p * i) + k]);
				}
			break;
			case 1:
				for (; i < m; i++) {
#	if MC_TARGET_HAVE_FMA
					s = mc_fmal(a[(n * i) + j], b[(p * i) + k], s);
#	else
//!# Accurate dot product sum(x[i] * y[i], i=0...n-1) of two vectors.
//!# Accurate Sum and Dot Product, Takeshi Ogita, Siegfried M. Rump
//!# and Shin'ichi Oishi 2005, published in SIAM Journal on Scientific
//!# Computing (SISC), 26(6):1955-1988, 2005.

//!# TwoProduct(a[(n * i) + j],b[(p * i) + k],h,r).
					q  = a[(n * i) + j];
//!# split a[(n * i) + j] into x1,x2.
					r  = cs * q;
					x2 = r - q;
					x1 = r - x2;
					x2 = q - x1;
					r  = b[(p * i) + k];
//!# h=a[(n * i) + j]*b[(p * i) + k].
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

#endif /* !MC_DOTPMX1_H */

/* EOF */