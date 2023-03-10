//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zgammaln.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_atan.h>
#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_fisint.h>
#include <macadam/details/math/mc_itrunc.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_pow.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_sinh.h>
#include <macadam/details/math/mc_zmod.h>
#include <macadam/details/math/mc_znorm.h>

#ifndef MC_ZLGAMMA_H
#define MC_ZLGAMMA_H

#pragma mark - mc_zgammaln_approx0 -

MC_TARGET_PROC void mc_zgammalnf_approx0(float * r_r, float * r_i, float x_r, float x_i)
{
//!# Computes the log of gamma function for complex argument.
//!#
//!# \note: Shanjie Zhang, Jianming Jin, Computation of Special
//!# Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45.
//!#
//!# This routine is copyrighted by Shanjie Zhang and Jianming Jin.
//!# However, they give permission to incorporate this routine into
//!# a user program provided that the copyright is acknowledged.
//!#
//!# \note: However, we don't know what it means as this is a full
//!# rewrite from Fortran to C; + extremely greedy, inelegant to
//!# the least. Todo: something usable.
//!#
	const float a1  = +0.08333333333333333000000000000000000000E+00f;
	const float a2  = -2.77777777777777800000000000000000000000E-03f;
	const float a3  = +7.93650793650793700000000000000000000000E-04f;
	const float a4  = -5.95238095238095200000000000000000000000E-04f;
	const float a5  = +8.41750841750841800000000000000000000000E-04f;
	const float a6  = -1.91752691752691800000000000000000000000E-03f;
	const float a7  = +6.41025641025641000000000000000000000000E-03f;
	const float a8  = -0.02955065359477124000000000000000000000E+00f;
	const float a9  = +0.17964437236883070000000000000000000000E+00f;
	const float a10 = -1.39243221690590000000000000000000000000E+00f;

	float a, b, c, d;
	int n = 0, i;

	*r_r = -MCLIMITS_MAXF;
	*r_i = +0.0f;

	if (!(x_i == 0.0f && mc_fisintf(x_r) && x_r <= 0.0f)) {
		c = x_r;
		if (x_r < 0.0f) {
			x_r = -x_r;
			x_i = -x_i;
		}
		d = x_r;
		if (x_r < 7.0f) {
			n = mc_itruncf(7.0f - x_r); 
			d = x_r + mc_cast(float, n);
		}
		 a   = mc_zmodf(d, x_i);
		 b   = mc_atanf(x_i / d);
		*r_r = (d - 0.5f) * mc_logf(a) - b * x_i - d + 0.5f * mc_logf(MCK_KF(MCK_2PI));
		*r_i = b * (d - 0.5f) + x_i * mc_logf(a) - x_i;

		 d   = mc_powf(a, -1.0f);
		*r_r = *r_r + (a1 * d * mc_cosf(b));
		*r_i = *r_i - (a1 * d * mc_sinf(b));

		 d   = mc_powf(a, -3.0f);
		*r_r = *r_r + (a2 * d * mc_cosf(3.0f * b));
		*r_i = *r_i - (a2 * d * mc_sinf(3.0f * b));

		 d   = mc_powf(a, -5.0f);
		*r_r = *r_r + (a3 * d * mc_cosf(5.0f * b));
		*r_i = *r_i - (a3 * d * mc_sinf(5.0f * b));

		 d   = mc_powf(a, -7.0f);
		*r_r = *r_r + (a4 * d * mc_cosf(7.0f * b));
		*r_i = *r_i - (a4 * d * mc_sinf(7.0f * b));

		 d   = mc_powf(a, -9.0f);
		*r_r = *r_r + (a5 * d * mc_cosf(9.0f * b));
		*r_i = *r_i - (a5 * d * mc_sinf(9.0f * b));

		 d   = mc_powf(a, -11.0f);
		*r_r = *r_r + (a6 * d * mc_cosf(11.0f * b));
		*r_i = *r_i - (a6 * d * mc_sinf(11.0f * b));

		 d   = mc_powf(a, -13.0f);
		*r_r = *r_r + (a7 * d * mc_cosf(13.0f * b));
		*r_i = *r_i - (a7 * d * mc_sinf(13.0f * b));

		 d   = mc_powf(a, -15.0f);
		*r_r = *r_r + (a8 * d * mc_cosf(15.0f * b));
		*r_i = *r_i - (a8 * d * mc_sinf(15.0f * b));

		 d   = mc_powf(a, -17.0f);
		*r_r = *r_r + (a9 * d * mc_cosf(17.0f * b));
		*r_i = *r_i - (a9 * d * mc_sinf(17.0f * b));

		 d   = mc_powf(a, -19.0f);
		*r_r = *r_r + (a10 * d * mc_cosf(19.0f * b));
		*r_i = *r_i - (a10 * d * mc_sinf(19.0f * b));

		if (x_r <= 7.0f) {
			b = 0.0f;
			a = 0.0f;
			for (i = 0 ; i <= (n - 1); i++) {
				d = x_r + mc_cast(float, i);
				b = b + (0.5f * mc_logf(mc_znormf(d, x_i))); 
				a = a + (mc_atanf(x_i / (x_r + mc_cast(float, i))));
			}
			*r_r = *r_r - b;
			*r_i = *r_i - a;
		}
		if (c < 0.0f) {
			 a   =  mc_zmodf(x_r, x_i);
			 c   = -mc_sinf(MCK_KF(MCK_PI) * x_r) * mc_coshf(MCK_KF(MCK_PI) * x_i);
			 d   = -mc_cosf(MCK_KF(MCK_PI) * x_r) * mc_sinhf(MCK_KF(MCK_PI) * x_i);
			 b   =  mc_atanf(d / c);
			 b   =  c < 0.0f ? b + MCK_KF(MCK_PI) : b;
			*r_r =  mc_logf(MCK_KF(MCK_PI) / (a * mc_zmodf(c, d))) - (*r_r);
			*r_i = -mc_atanf(x_i / x_r) - b - (*r_i);
		}
	}
}

MC_TARGET_PROC void mc_zgammaln_approx0(double * r_r, double * r_i, double x_r, double x_i)
{
//!# Computes the log of gamma function for complex argument.
//!#
//!# \note: Shanjie Zhang, Jianming Jin, Computation of Special
//!# Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45.
//!#
//!# This routine is copyrighted by Shanjie Zhang and Jianming Jin.
//!# However, they give permission to incorporate this routine into
//!# a user program provided that the copyright is acknowledged.
//!#
//!# \note: However, we don't know what it means as this is a full
//!# rewrite from Fortran to C; + extremely greedy, inelegant to
//!# the least. Todo: something usable.
//!#
	const double a1  = +0.0833333333333333300000000000000000000000E+00;
	const double a2  = -2.7777777777777780000000000000000000000000E-03;
	const double a3  = +7.9365079365079370000000000000000000000000E-04;
	const double a4  = -5.9523809523809520000000000000000000000000E-04;
	const double a5  = +8.4175084175084180000000000000000000000000E-04;
	const double a6  = -1.9175269175269180000000000000000000000000E-03;
	const double a7  = +6.4102564102564100000000000000000000000000E-03;
	const double a8  = -0.0295506535947712400000000000000000000000E+00;
	const double a9  = +0.1796443723688307000000000000000000000000E+00;
	const double a10 = -1.3924322169059000000000000000000000000000E+00;

	double a, b, c, d;
	int n = 0, i;

	*r_r = -MCLIMITS_MAX;
	*r_i = +0.0;

	if (!(x_i == 0.0 && mc_fisint(x_r) && x_r <= 0.0)) {
		c = x_r;
		if (x_r < 0.0) {
			x_r = -x_r;
			x_i = -x_i;
		}
		d = x_r;
		if (x_r < 7.0) {
			n = mc_itrunc(7.0 - x_r); 
			d = x_r + mc_cast(double, n);
		}
		 a   = mc_zmod(d, x_i);
		 b   = mc_atan(x_i / d);
		*r_r = (d - 0.5) * mc_log(a) - b * x_i - d + 0.5 * mc_log(MCK_K(MCK_2PI));
		*r_i = b * (d - 0.5) + x_i * mc_log(a) - x_i;

		 d   = mc_pow(a, -1.0);
		*r_r = *r_r + (a1 * d * mc_cos(b));
		*r_i = *r_i - (a1 * d * mc_sin(b));

		 d   = mc_pow(a, -3.0);
		*r_r = *r_r + (a2 * d * mc_cos(3.0 * b));
		*r_i = *r_i - (a2 * d * mc_sin(3.0 * b));

		 d   = mc_pow(a, -5.0);
		*r_r = *r_r + (a3 * d * mc_cos(5.0 * b));
		*r_i = *r_i - (a3 * d * mc_sin(5.0 * b));

		 d   = mc_pow(a, -7.0);
		*r_r = *r_r + (a4 * d * mc_cos(7.0 * b));
		*r_i = *r_i - (a4 * d * mc_sin(7.0 * b));

		 d   = mc_pow(a, -9.0);
		*r_r = *r_r + (a5 * d * mc_cos(9.0 * b));
		*r_i = *r_i - (a5 * d * mc_sin(9.0 * b));

		 d   = mc_pow(a, -11.0);
		*r_r = *r_r + (a6 * d * mc_cos(11.0 * b));
		*r_i = *r_i - (a6 * d * mc_sin(11.0 * b));

		 d   = mc_pow(a, -13.0);
		*r_r = *r_r + (a7 * d * mc_cos(13.0 * b));
		*r_i = *r_i - (a7 * d * mc_sin(13.0 * b));

		 d   = mc_pow(a, -15.0);
		*r_r = *r_r + (a8 * d * mc_cos(15.0 * b));
		*r_i = *r_i - (a8 * d * mc_sin(15.0 * b));

		 d   = mc_pow(a, -17.0);
		*r_r = *r_r + (a9 * d * mc_cos(17.0 * b));
		*r_i = *r_i - (a9 * d * mc_sin(17.0 * b));

		 d   = mc_pow(a, -19.0);
		*r_r = *r_r + (a10 * d * mc_cos(19.0 * b));
		*r_i = *r_i - (a10 * d * mc_sin(19.0 * b));

		if (x_r <= 7.0) {
			b = 0.0;
			a = 0.0;
			for (i = 0; i <= (n - 1); i++) {
				d = x_r + mc_cast(double, i);
				b = b + (0.5 * mc_log(mc_znorm(d, x_i))); 
				a = a + (mc_atan(x_i / (x_r + mc_cast(double, i))));
			} 
			*r_r = *r_r - b;
			*r_i = *r_i - a;
		}
		if (c < 0.0) {
			 a   =  mc_zmod(x_r, x_i);
			 c   = -mc_sin(MCK_K(MCK_PI) * x_r) * mc_cosh(MCK_K(MCK_PI) * x_i);
			 d   = -mc_cos(MCK_K(MCK_PI) * x_r) * mc_sinh(MCK_K(MCK_PI) * x_i);
			 b   =  mc_atan(d / c);
			 b   =  c < 0.0 ? b + MCK_K(MCK_PI) : b;
			*r_r =  mc_log(MCK_K(MCK_PI) / (a *mc_zmod(c, d))) - (*r_r);
			*r_i = -mc_atan(x_i / x_r) - b - (*r_i);
		}
	}
}

MC_TARGET_PROC void mc_zgammalnl_approx0(long double * r_r, long double * r_i, long double x_r, long double x_i)
{
//!# Computes the log of gamma function for complex argument.
//!#
//!# \note: Shanjie Zhang, Jianming Jin, Computation of Special
//!# Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45.
//!#
//!# This routine is copyrighted by Shanjie Zhang and Jianming Jin.
//!# However, they give permission to incorporate this routine into
//!# a user program provided that the copyright is acknowledged.
//!#
//!# \note: However, we don't know what it means as this is a full
//!# rewrite from Fortran to C; + extremely greedy, inelegant to
//!# the least. Todo: something usable.
//!#
	const long double a1  = +0.083333333333333330000000000000000000000000000000000000000000000E+00L;
	const long double a2  = -2.777777777777778000000000000000000000000000000000000000000000000E-03L;
	const long double a3  = +7.936507936507937000000000000000000000000000000000000000000000000E-04L;
	const long double a4  = -5.952380952380952000000000000000000000000000000000000000000000000E-04L;
	const long double a5  = +8.417508417508418000000000000000000000000000000000000000000000000E-04L;
	const long double a6  = -1.917526917526918000000000000000000000000000000000000000000000000E-03L;
	const long double a7  = +6.410256410256410000000000000000000000000000000000000000000000000E-03L;
	const long double a8  = -0.029550653594771240000000000000000000000000000000000000000000000E+00L;
	const long double a9  = +0.179644372368830700000000000000000000000000000000000000000000000E+00L;
	const long double a10 = -1.392432216905900000000000000000000000000000000000000000000000000E+00L;

	long double a, b, c, d;
	int n = 0, i;

	*r_r = -MCLIMITS_MAXL;
	*r_i = +0.0L;

	if (!(x_i == 0.0L && mc_fisintl(x_r) && x_r <= 0.0L)) {
		c = x_r;
		if (x_r < 0.0L) {
			x_r = -x_r;
			x_i = -x_i;
		}
		d = x_r;
		if (x_r < 7.0L) {
			n = mc_itruncl(7.0L - x_r); 
			d = x_r + mc_cast(long double, n);
		}
		 a   = mc_zmodl(d, x_i);
		 b   = mc_atanl(x_i / d);
		*r_r = (d - 0.5L) * mc_logl(a) - b * x_i - d + 0.5L * mc_logl(MCK_KL(MCK_2PI));
		*r_i = b * (d - 0.5L) + x_i * mc_logl(a) - x_i;

		 d   = mc_powl(a, -1.0L);
		*r_r = *r_r + (a1 * d * mc_cosl(b));
		*r_i = *r_i - (a1 * d * mc_sinl(b));

		 d   = mc_powl(a, -3.0L);
		*r_r = *r_r + (a2 * d * mc_cosl(3.0L * b));
		*r_i = *r_i - (a2 * d * mc_sinl(3.0L * b));

		 d   = mc_powl(a, -5.0L);
		*r_r = *r_r + (a3 * d * mc_cosl(5.0L * b));
		*r_i = *r_i - (a3 * d * mc_sinl(5.0L * b));

		 d   = mc_powl(a, -7.0L);
		*r_r = *r_r + (a4 * d * mc_cosl(7.0L * b));
		*r_i = *r_i - (a4 * d * mc_sinl(7.0L * b));

		 d   = mc_powl(a, -9.0L);
		*r_r = *r_r + (a5 * d * mc_cosl(9.0L * b));
		*r_i = *r_i - (a5 * d * mc_sinl(9.0L * b));

		 d   = mc_powl(a, -11.0L);
		*r_r = *r_r + (a6 * d * mc_cosl(11.0L * b));
		*r_i = *r_i - (a6 * d * mc_sinl(11.0L * b));

		 d   = mc_powl(a, -13.0L);
		*r_r = *r_r + (a7 * d * mc_cosl(13.0L * b));
		*r_i = *r_i - (a7 * d * mc_sinl(13.0L * b));

		 d   = mc_powl(a, -15.0L);
		*r_r = *r_r + (a8 * d * mc_cosl(15.0L * b));
		*r_i = *r_i - (a8 * d * mc_sinl(15.0L * b));

		 d   = mc_powl(a, -17.0L);
		*r_r = *r_r + (a9 * d * mc_cosl(17.0L * b));
		*r_i = *r_i - (a9 * d * mc_sinl(17.0L * b));

		 d   = mc_powl(a, -19.0L);
		*r_r = *r_r + (a10 * d * mc_cosl(19.0L * b));
		*r_i = *r_i - (a10 * d * mc_sinl(19.0L * b));

		if (x_r <= 7.0L) {
			b = 0.0L;
			a = 0.0L;
			for (i = 0; i <= (n - 1); i++) {
				d = x_r + mc_cast(long double, i);
				b = b + (0.5L * mc_logl(mc_znorml(d, x_i))); 
				a = a + (mc_atanl(x_i / (x_r + mc_cast(long double, i))));
			} 
			*r_r = *r_r - b;
			*r_i = *r_i - a;
		}
		if (c < 0.0L) {
			a    =  mc_zmodl(x_r, x_i);
			c    = -mc_sinl(MCK_KL(MCK_PI) * x_r) * mc_coshl(MCK_KL(MCK_PI) * x_i);
			d    = -mc_cosl(MCK_KL(MCK_PI) * x_r) * mc_sinhl(MCK_KL(MCK_PI) * x_i);
			b    =  mc_atanl(d / c);
			b    =  c < 0.0L ? b + MCK_KL(MCK_PI) : b;
			*r_r =  mc_logl(MCK_KL(MCK_PI) / (a * mc_zmodl(c, d))) - (*r_r);
			*r_i = -mc_atanl(x_i / x_r) - b - (*r_i);
		}
	}
}

#pragma mark - mc_zgammaln -

MC_TARGET_PROC void mc_zgammalnf(float * r_r, float * r_i, const float x_r, const float x_i)
{
	mc_zgammalnf_approx0(r_r, r_i, x_r, x_i);
}

MC_TARGET_PROC void mc_zgammaln(double * r_r, double * r_i, const double x_r, const double x_i)
{
	mc_zgammaln_approx0(r_r, r_i, x_r, x_i);
}

MC_TARGET_PROC void mc_zgammalnl(long double * r_r, long double * r_i, const long double x_r, const long double x_i)
{
	mc_zgammalnl_approx0(r_r, r_i, x_r, x_i);
}

#endif /* !MC_ZLGAMMA_H */

/* EOF */