//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_igamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_erf.h>
#include <macadam/details/math/mc_erfc.h>
#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_fisint.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_factorial.h>
#include <macadam/details/math/mc_gamma.h>
#include <macadam/details/math/mc_gammaln.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_pow.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_IGAMMA_H
#define MC_IGAMMA_H

#pragma mark - mc_igamma_pseries_approx0 -

MC_TARGET_PROC float mc_igamma_pseriesf_approx0(const float a, const float z)
{
//!# Power series.
	const float e = MCLIMITS_EPSILONF;
	float term    = 1.0f / a;
	float sum     = term;
	float w       = a;
	do {
		w    = w + 1.0f;
		term = term * (z / w);
		sum  = sum + term;
	} while (mc_fabsf(term) > mc_fabsf(sum) * e);
	return sum;
}

MC_TARGET_PROC double mc_igamma_pseries_approx0(const double a, const double z)
{
//!# Power series.
	const double e = MCLIMITS_EPSILON;
	double term    = 1.0 / a;
	double sum     = term;
	double w       = a;
	do {
		w    = w + 1.0;
		term = term * (z / w);
		sum  = sum + term;
	} while (mc_fabs(term) > mc_fabs(sum) * e);
	return sum;
}

MC_TARGET_PROC long double mc_igamma_pseriesl_approx0(const long double a, const long double z)
{
//!# Power series.
	const long double e = MCLIMITS_EPSILONL;
	long double term    = 1.0L / a;
	long double sum     = term;
	long double w       = a;
	do {
		w    = w + 1.0L;
		term = term * (z / w);
		sum  = sum + term;
	} while (mc_fabsl(term) > mc_fabsl(sum) * e);
	return sum;
}

#pragma mark - mc_igamma_cfrac_approx0 -

MC_TARGET_PROC float mc_igamma_cfracf_approx0(const float a, const float z)
{
//!# Continued fraction.
	const float e1 = MCLIMITS_EPSILONF;
	const float e3 = e1 * 3.0f;
	float k = 1.0f, c0, c1, c2, c3, c4, c5;
	c1 = z + 1.0f - a;
	c2 = 1.0f / e3;
	c3 = 1.0f / c1;
	c5 = c3;
	do {
		c0 = k * (a - k);
		c1 = c1 +  2.0f;
		c3 = c0 * c3 + c1;
		c3 = mc_fabsf(c3) < e3 ? e3 : c3;
		c2 = c1 + c0 / c2;
		c2 = mc_fabsf(c2) < e3 ? e3 : c2;
		c3 = 1.0f / c3;
		c4 = c3 * c2;
		c5 = c5 * c4;
		k  = k + 1.0f;
	} while (mc_fabsf(c4 - 1.0f) > e1);
	return c5;
}

MC_TARGET_PROC double mc_igamma_cfrac_approx0(const double a, const double z)
{
//!# Continued fraction.
	const double e1 = MCLIMITS_EPSILON;
	const double e3 = e1 * 3.0;
	double k = 1.0, c0, c1, c2, c3, c4, c5;
	c1 = z + 1.0 - a;
	c2 = 1.0 / e3;
	c3 = 1.0 / c1;
	c5 = c3;
	do {
		c0 = k * (a - k);
		c1 = c1 +  2.0;
		c3 = c0 * c3 + c1;
		c3 = mc_fabs(c3) < e3 ? e3 : c3;
		c2 = c1 + c0 / c2;
		c2 = mc_fabs(c2) < e3 ? e3 : c2;
		c3 = 1.0 / c3;
		c4 = c3 * c2;
		c5 = c5 * c4;
		k   = k + 1.0;
	} while (mc_fabs(c4 - 1.0) > e1);
	return c5;
}

MC_TARGET_PROC long double mc_igamma_cfracl_approx0(const long double a, const long double z)
{
//!# Continued fraction.
	const long double e1 = MCLIMITS_EPSILONL;
	const long double e3 = e1 * 3.0L;
	long double k = 1.0L, c0, c1, c2, c3, c4, c5;
	c1 = z + 1.0L - a;
	c2 = 1.0L / e3;
	c3 = 1.0L / c1;
	c5 = c3;
	do {
		c0 = k * (a - k);
		c1 = c1 +  2.0L;
		c3 = c0 * c3 + c1;
		c3 = mc_fabsl(c3) < e3 ? e3 : c3;
		c2 = c1 + c0 / c2;
		c2 = mc_fabsl(c2) < e3 ? e3 : c2;
		c3 = 1.0L / c3;
		c4 = c3 * c2;
		c5 = c5 * c4;
		k  = k + 1.0L;
	} while (mc_fabsl(c4 - 1.0L) > e1);
	return c5;
}

#pragma mark - mc_igamma_lower -

MC_TARGET_PROC float mc_igamma_lowerf_approx0(const float a, const float z)
{
//!# Lower incomplete gamma function of a and z.
	if (a > 0.0f && z > 0.0f) {
		const float y = a * mc_logf(z) - z;
		if (y < -FLT_MAX_10_EXP) {
			return 0.0f;
		}
		return mc_igamma_pseriesf_approx0(a, z) * mc_expf(y);
	}
	return MCK_NAN;
}

MC_TARGET_PROC double mc_igamma_lower_approx0(const double a, const double z)
{
//!# Lower incomplete gamma function of a and z.
	if (a > 0.0 && z > 0.0) {
		const double y = a * mc_log(z) - z;
		if (y < -DBL_MAX_10_EXP) {
			return 0.0;
		}
		return mc_igamma_pseries_approx0(a, z) * mc_exp(y);
	}
	return MCK_NAN;
}

MC_TARGET_PROC long double mc_igamma_lowerl_approx0(const long double a, const long double z)
{
//!# Lower incomplete gamma function of a and z.
	if (a > 0.0L && z > 0.0L) {
		const long double y = a * mc_logl(z) - z;
		if (y < -LDBL_MAX_10_EXP) {
			return 0.0L;
		}
		return mc_igamma_pseriesl_approx0(a, z) * mc_expl(y);
	}
	return MCK_NAN;
}

#pragma mark - mc_igamma_upper_approx0 -

MC_TARGET_PROC float mc_igamma_upperf_approx0(const float a, const float z)
{
//!# Upper incomplete gamma function of a and z.
	if (a > 0.0f && z > 0.0f) {
		const float w = a * mc_logf(z) - z;
		if (w < -FLT_MAX_10_EXP) {
			return 0.0f;
		}
		return mc_igamma_cfracf_approx0(a, z) * mc_expf(w);
	}
	return MCK_NAN;
}

MC_TARGET_PROC double mc_igamma_upper_approx0(const double a, const double z)
{
//!# Upper incomplete gamma function of a and z.
	if (a > 0.0 && z > 0.0) {
		const double w = a * mc_log(z) - z;
		if (w < -DBL_MAX_10_EXP) {
			return 0.0;
		}
		return mc_igamma_cfrac_approx0(a, z) * mc_exp(w);
	}
	return MCK_NAN;
}

MC_TARGET_PROC long double mc_igamma_upperl_approx0(const long double a, const long double z)
{
//!# Upper incomplete gamma function of a and z.
	if (a > 0.0L && z > 0.0L) {
		const long double w = a * mc_logl(z) - z;
		if (w < -LDBL_MAX_10_EXP) {
			return 0.0L;
		}
		return mc_igamma_cfracl_approx0(a, z) * mc_expl(w);
	}
	return MCK_NAN;
}

#pragma mark - mc_igamma_p_approx0 -

MC_TARGET_PROC float mc_igamma_pf_approx0(const float a, const float z)
{
//!# Normalised lower incomplete gamma function of a and z.
	float p = MCK_NAN;
	if (a > 0.0f && z >= 0.0f) {
		if (z <= 0.0f) {
			return 0.0f;
		}
		const float y = a * mc_logf(z) - z;
		if (y >= -FLT_MAX_10_EXP) {
			const float w = y - mc_gammalnf_approx0(a, MC_NULLPTR);
			if (z < a + 1.0f) {
				p = mc_igamma_pseriesf_approx0(a, z) * mc_expf(w);
			} else {
				p = 1.0f - mc_igamma_cfracf_approx0(a, z) * mc_expf(w);
			}
		} else {
			p = 0.0f;
		}
	}
	return p;
}

MC_TARGET_PROC double mc_igamma_p_approx0(const double a, const double z)
{
//!# Normalised lower incomplete gamma function of a and z.
	double p = MCK_NAN;
	if (a > 0.0 && z >= 0.0) {
		if (z <= 0.0) {
			return 0.0;
		}
		const double y = a * mc_log(z) - z;
		if (y >= -DBL_MAX_10_EXP) {
			const double w = y - mc_gammaln_approx0(a, MC_NULLPTR);
			if (z < a + 1.0) {
				p = mc_igamma_pseries_approx0(a, z) * mc_exp(w);
			} else {
				p = 1.0 - mc_igamma_cfrac_approx0(a, z) * mc_exp(w);
			}
		} else {
			p = 0.0;
		}
	}
	return p;
}

MC_TARGET_PROC long double mc_igamma_pl_approx0(const long double a, const long double z)
{
//!# Normalised lower incomplete gamma function of a and z.
	long double p = MCK_NAN;
	if (a > 0.0L && z >= 0.0L) {
		if (z <= 0.0L) {
			return 0.0L;
		}
		const long double y = a * mc_logl(z) - z;
		if (y >= -LDBL_MAX_10_EXP) {
			const long double w = y - mc_gammalnl_approx0(a, MC_NULLPTR);
			if (z < a + 1.0L) {
				p = mc_igamma_pseriesl_approx0(a, z) * mc_expl(w);
			} else {
				p = 1.0L - mc_igamma_cfracl_approx0(a, z) * mc_expl(w);
			}
		} else {
			p = 0.0L;
		}
	}
	return p;
}

#pragma mark - mc_igamma_q_approx0 -

MC_TARGET_PROC float mc_igamma_qf_approx0(const float a, const float z)
{
//!# Normalised upper incomplete gamma function of a and z.
	float q = MCK_NAN;
	if (a > 0.0f && z >= 0.0f) {
		if (z <= 0.0f) {
			return 1.0f;
		}
		const float y = a * mc_logf(z) - z;
		if (y >= -FLT_MAX_10_EXP) {
			const float w = y - mc_gammalnf_approx0(a, MC_NULLPTR);
			if (z >= a + 1.0f) {
				q = mc_igamma_cfracf_approx0(a, z) * mc_expf(w);
			} else {
				q = 1.0f - mc_igamma_pseriesf_approx0(a, z) * mc_expf(w);
			}
		} else {
			q = 1.0f;
		}
	}
	return q;
}

MC_TARGET_PROC double mc_igamma_q_approx0(const double a, const double z)
{
//!# Normalised upper incomplete gamma function of a and z.
	double q = MCK_NAN;
	if (a > 0.0 && z >= 0.0) {
		if (z <= 0.0) {
			return 1.0;
		}
		const double y = a * mc_log(z) - z;
		if (y >= -DBL_MAX_10_EXP) {
			const double w = y - mc_gammaln_approx0(a, MC_NULLPTR);
			if (z >= a + 1.0) {
				q = mc_igamma_cfrac_approx0(a, z) * mc_exp(w);
			} else {
				q = 1.0 - mc_igamma_pseries_approx0(a, z) * mc_exp(w);
			}
		} else {
			q = 1.0;
		}
	}
	return q;
}

MC_TARGET_PROC long double mc_igamma_ql_approx0(const long double a, const long double z)
{
//!# Normalised upper incomplete gamma function of a and z.
	long double q = MCK_NAN;
	if (a > 0.0L && z >= 0.0L) {
		if (z <= 0.0L) {
			return 1.0L;
		}
		const long double y = a * mc_logl(z) - z;
		if (y >= -LDBL_MAX_10_EXP) {
			const long double w = y - mc_gammalnl_approx0(a, MC_NULLPTR);
			if (z >= a + 1.0L) {
				q = mc_igamma_cfracl_approx0(a, z) * mc_expl(w);
			} else {
				q = 1.0L - mc_igamma_pseriesl_approx0(a, z) * mc_expl(w);
			}
		} else {
			q = 1.0L;
		}
	}
	return q;
}

#pragma mark - mc_igamma_small_approx0 -

MC_TARGET_PROC float mc_igamma_smallf_approx0(const float a, const float z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| <= 1.
	float r = 0.0f, p;
	if (a > 0.0f) {
		p = 1.0f;

		r = r + (p / mc_factorialf(0.0f))  / (0.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(1.0f))  / (1.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(2.0f))  / (2.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(3.0f))  / (3.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(4.0f))  / (4.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(5.0f))  / (5.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(6.0f))  / (6.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(7.0f))  / (7.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(8.0f))  / (8.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(9.0f))  / (9.0f  + a);
		p = p * -z;
		r = r + (p / mc_factorialf(10.0f)) / (10.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(11.0f)) / (11.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(12.0f)) / (12.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(13.0f)) / (13.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(14.0f)) / (14.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(15.0f)) / (15.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(16.0f)) / (16.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(17.0f)) / (17.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(18.0f)) / (18.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(19.0f)) / (19.0f + a);
		p = p * -z;
		r = r + (p / mc_factorialf(20.0f)) / (20.0f + a);

		r = mc_expf(a * mc_logf(z) + mc_logf(r) - mc_gammalnf_approx0(a, MC_NULLPTR));
	}
	return r;
}

MC_TARGET_PROC double mc_igamma_small_approx0(const double a, const double z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| <= 1.
	double r = 0.0, p;
	if (a > 0.0) {
		p = 1.0;

		r = r + (p / mc_factorial(0.0))  / (0.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(1.0))  / (1.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(2.0))  / (2.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(3.0))  / (3.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(4.0))  / (4.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(5.0))  / (5.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(6.0))  / (6.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(7.0))  / (7.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(8.0))  / (8.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(9.0))  / (9.0  + a);
		p = p * -z;
		r = r + (p / mc_factorial(10.0)) / (10.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(11.0)) / (11.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(12.0)) / (12.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(13.0)) / (13.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(14.0)) / (14.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(15.0)) / (15.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(16.0)) / (16.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(17.0)) / (17.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(18.0)) / (18.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(19.0)) / (19.0 + a);
		p = p * -z;
		r = r + (p / mc_factorial(20.0)) / (20.0 + a);

		r = mc_exp(a * mc_log(z) + mc_log(r) - mc_gammaln_approx0(a, MC_NULLPTR));
	}
	return r;
}

MC_TARGET_PROC long double mc_igamma_smalll_approx0(const long double a, const long double z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| <= 1.
	long double r = 0.0L, p;
	if (a > 0.0L) {
		p = 1.0L;

		r = r + (p / mc_factoriall(0.0L))  / (0.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(1.0L))  / (1.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(2.0L))  / (2.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(3.0L))  / (3.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(4.0L))  / (4.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(5.0L))  / (5.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(6.0L))  / (6.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(7.0L))  / (7.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(8.0L))  / (8.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(9.0L))  / (9.0L  + a);
		p = p * -z;
		r = r + (p / mc_factoriall(10.0L)) / (10.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(11.0L)) / (11.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(12.0L)) / (12.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(13.0L)) / (13.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(14.0L)) / (14.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(15.0L)) / (15.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(16.0L)) / (16.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(17.0L)) / (17.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(18.0L)) / (18.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(19.0L)) / (19.0L + a);
		p = p * -z;
		r = r + (p / mc_factoriall(20.0L)) / (20.0L + a);

		r = mc_expl(a * mc_logl(z) + mc_logl(r) - mc_gammalnl_approx0(a, MC_NULLPTR));
	}
	return r;
}

#pragma mark - mc_igamma_medium_approx0 -

MC_TARGET_PROC float mc_igamma_mediumf_approx0(const float a, const float z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| < a + 1.
	float r = 0.0f;
	float p = 1.0f / a;
	float c = p;
	float s = p;
	float d = -s + c;
	float l = s;
	float m = 0.0f;
	float e = 0.0f;
	int k   = 1;
	if (a > 0.0f) {
		const float y = a * mc_logf(z) - z;
		if (y < -FLT_MAX_10_EXP) {
			return r;
		}
		r = mc_expf(y - mc_gammalnf_approx0(a, MC_NULLPTR));
		if (r > 0.0f) {
			e = 1.0E-10f / r;
		}
		if (e <= 0.0f) {
			e = 1.0E-10f;
		}
		for (; p > e * l; k++) {
			p = p * ((1.0f / z) * (a + mc_cast(float, k)));
			c = p + d;
			s = l + c;
			d = (l - s) + c;
			l = s;
		}
		m = l;
		l = l * r;
		d = d + (1.0f / (m - l) * r);
		p = p * ((1.0f / z) * (a + mc_cast(float, k)));
		m = p + d;
		k = k + 1;
		for (; (p + d) > e * m; k++) {
			p = p * ((1.0f / z) * (a + mc_cast(float, k)));
			c = p + d;
			s = m + c;
			d = (m - s) + c;
			m = s;
		}
		r = l + ((m + d) * r);
	}
	return r;
}

MC_TARGET_PROC double mc_igamma_medium_approx0(const double a, const double z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| < a + 1.
	double r = 0.0;
	double p = 1.0 / a;
	double c = p;
	double s = p;
	double d = -s + c;
	double l = s;
	double m = 0.0;
	double e = 0.0;
	int k   = 1;
	if (a > 0.0) {
		const double y = a * mc_log(z) - z;
		if (y < -DBL_MAX_10_EXP) {
			return r;
		}
		r = mc_exp(y - mc_gammaln_approx0(a, MC_NULLPTR));
		if (r > 0.0) {
			e = 1.0E-15 / r;
		}
		if (e <= 0.0) {
			e = 1.0E-15;
		}
		for (; p > e * l; k++) {
			p = p * ((1.0 / z) * (a + mc_cast(double, k)));
			c = p + d;
			s = l + c;
			d = (l - s) + c;
			l = s;
		}
		m = l;
		l = l * r;
		d = d + (1.0 / (m - l) * r);
		p = p * ((1.0 / z) * (a + mc_cast(double, k)));
		m = p + d;
		k = k + 1;
		for (; (p + d) > e * m; k++) {
			p = p * ((1.0 / z) * (a + mc_cast(double, k)));
			c = p + d;
			s = m + c;
			d = (m - s) + c;
			m = s;
		}
		r = l + ((m + d) * r);
	}
	return r;
}

MC_TARGET_PROC long double mc_igamma_mediuml_approx0(const long double a, const long double z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| < a + 1.
	long double r = 0.0L;
	long double p = 1.0L / a;
	long double c = p;
	long double s = p;
	long double d = -s + c;
	long double l = s;
	long double m = 0.0L;
	long double e = 0.0L;
	int k   = 1;
	if (a > 0.0L) {
		const long double y = a * mc_logl(z) - z;
		if (y < -LDBL_MAX_10_EXP) {
			return r;
		}
		r = mc_expl(y - mc_gammalnl_approx0(a, MC_NULLPTR));
		if (r > 0.0L) {
			e = 1.0E-15L / r;
		}
		if (e <= 0.0L) {
			e = 1.0E-15L;
		}
		for (; p > e * l; k++) {
			p = p * ((1.0L / z) * (a + mc_cast(long double, k)));
			c = p + d;
			s = l + c;
			d = (l - s) + c;
			l = s;
		}
		m = l;
		l = l * r;
		d = d + (1.0L / (m - l) * r);
		p = p * ((1.0L / z) * (a + mc_cast(long double, k)));
		m = p + d;
		k = k + 1;
		for (; (p + d) > e * m; k++) {
			p = p * ((1.0L / z) * (a + mc_cast(long double, k)));
			c = p + d;
			s = m + c;
			d = (m - s) + c;
			m = s;
		}
		r = l + ((m + d) * r);
	}
	return r;
}

#pragma mark - mc_igamma_large_approx0 -

MC_TARGET_PROC float mc_igamma_largef_approx0(const float a, const float z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| >= a + 1.
	float r = 0.0f, p = 1.0f / a, s = p;
	int k = 1, n = mc_cast(int, (z - a - 1.0f)) + 1;
	if (a > 0.0f && n < MCLIMITS_IMAX) {
		for (; k < n; k++) {
			p = p * (1.0f / (z + mc_cast(float, k)));
			s = s + p;
		}
		r = mc_expf(mc_logf(s) + a * mc_logf(z) - z - mc_gammalnf_approx0(a, MC_NULLPTR)) + mc_igamma_smallf_approx0(a + mc_cast(float, n), z);
	}
	return r;
}

MC_TARGET_PROC double mc_igamma_large_approx0(const double a, const double z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| >= a + 1.
	double r = 0.0, p = 1.0 / a, s = p;
	int k = 1, n = mc_cast(int, (z - a - 1.0)) + 1;
	if (a > 0.0 && n < MCLIMITS_IMAX) {
		for (; k < n; k++) {
			p = p * (1.0 / (z + mc_cast(double, k)));
			s = s + p;
		}
		r = mc_exp(mc_log(s) + a * mc_log(z) - z - mc_gammaln_approx0(a, MC_NULLPTR)) + mc_igamma_small_approx0(a + mc_cast(double, n), z);
	}
	return r;
}

MC_TARGET_PROC long double mc_igamma_largel_approx0(const long double a, const long double z)
{
//!# Regularized incomplete gamma function; integral 0 to z with |z| >= a + 1.
	long double r = 0.0L, p = 1.0L / a, s = p;
	int k = 1, n = mc_cast(int, (z - a - 1.0L)) + 1;
	if (a > 0.0L && n < MCLIMITS_IMAX) {
		for (; k < n; k++) {
			p = p * (1.0L / (z + mc_cast(long double, k)));
			s = s + p;
		}
		r = mc_expl(mc_logl(s) + a * mc_logl(z) - z - mc_gammalnl_approx0(a, MC_NULLPTR)) + mc_igamma_smalll_approx0(a + mc_cast(long double, n), z);
	}
	return r;
}

#pragma mark - mc_igammae_approx0 -

MC_TARGET_PROC float mc_igammaef_approx0(const float z, const float a)
{
//!# Regularized incomplete gamma function; integral 0 to z.
	if (z == 0.0f || !(a > 0.0f)) {
		return 0.0f;
	}
	if (mc_fabsf(z) <= 1.0f) {
		return mc_igamma_smallf_approx0(a, z);
	}
	if (mc_fabsf(z) <= (a + 1.0f)) {
		return mc_igamma_mediumf_approx0(a, z);
	}
	return mc_igamma_largef_approx0(a, z);
}

MC_TARGET_PROC double mc_igammae_approx0(const double z, const double a)
{
//!# Regularized incomplete gamma function; integral 0 to z.
	if (z == 0.0 || !(a > 0.0)) {
		return 0.0;
	}
	if (mc_fabs(z) <= 1.0) {
		return mc_igamma_small_approx0(a, z);
	}
	if (mc_fabs(z) <= (a + 1.0)) {
		return mc_igamma_medium_approx0(a, z);
	}
	return mc_igamma_large_approx0(a, z);
}

MC_TARGET_PROC long double mc_igammael_approx0(const long double z, const long double a)
{
//!# Regularized incomplete gamma function; integral 0 to z.
	if (z == 0.0L || !(a > 0.0L)) {
		return 0.0L;
	}
	if (mc_fabsl(z) <= 1.0L) {
		return mc_igamma_smalll_approx0(a, z);
	}
	if (mc_fabsl(z) <= (a + 1.0L)) {
		return mc_igamma_mediuml_approx0(a, z);
	}
	return mc_igamma_largel_approx0(a, z);
}

#pragma mark - mc_igamma_pseries -

MC_TARGET_PROC float mc_igamma_pseriesf_approx1(const float a, const float z)
{
//!# Power series.
	const float e = MCLIMITS_EPSILONF;
	float term    = 1.0f / a;
	float sum     = term;
	float w       = a;
	for (;;) {
		w    = w + 1.0f;
		term = term * (z / w);
		sum  = sum + term;
		if (mc_fabsf(term) < mc_fabsf(sum) * e) {
			break;
		}
	}
	return sum;
}

MC_TARGET_PROC double mc_igamma_pseries_approx1(const double a, const double z)
{
//!# Power series.
	const double e = MCLIMITS_EPSILON;
	double term    = 1.0 / a;
	double sum     = term;
	double w       = a;
	for (;;) {
		w    = w + 1.0;
		term = term * (z / w);
		sum  = sum + term;
		if (mc_fabs(term) < mc_fabs(sum) * e) {
			break;
		}
	}
	return sum;
}

MC_TARGET_PROC long double mc_igamma_pseriesl_approx1(const long double a, const long double z)
{
//!# Power series.
	const long double e = MCLIMITS_EPSILONL;
	long double term    = 1.0 / a;
	long double sum     = term;
	long double w       = a;
	for (;;) {
		w    = w + 1.0L;
		term = term * (z / w);
		sum  = sum + term;
		if (mc_fabsl(term) < mc_fabsl(sum) * e) {
			break;
		}
	}
	return sum;
}

#pragma mark - mc_igamma_cfrac -

MC_TARGET_PROC float mc_igamma_cfracf_approx1(const float a, const float z)
{
	const float e1 = MCLIMITS_EPSILONF;
	const float e2 = e1 * 2.0f;
	float fa       = 1.0f - a + z;
	if (fa == 0.0f) {
		fa = e2;
	}
	float aa, ba, ca, c0 = fa, da, d0 = 0.0f, de = 0.0f;
	float k = 1.0f;
	do {
		ba = (1.0f + k * 2.0f) - a + z;
		aa = k * (a - k);
		da = ba + aa * d0;
		if (da == 0.0f) {
			da = e2;
		}
		ba = (1.0f + k * 2.0f) - a + z;
		aa = k * (a - k);
		ca = ba + aa / c0;
		if (ca == 0.0f) {
			ca = e2;
		}
		da = 1.0f / da;
		de = ca * da;
		fa = fa * de;
		d0 = da;
		c0 = ca;
		k  = k + 1.0f;
	} while (mc_fabsf(de - 1.0f) > e1);
	return fa;
}

MC_TARGET_PROC double mc_igamma_cfrac_approx1(const double a, const double z)
{
	const double e1 = MCLIMITS_EPSILON;
	const double e2 = e1 * 2.0;
	double fa       = 1.0 - a + z;
	if (fa == 0.0) {
		fa = e2;
	}
	double aa, ba, ca, c0 = fa, da, d0 = 0.0, de = 0.0;
	double k = 1.0;
	do {
		ba = (1.0 + k * 2.0) - a + z;
		aa = k * (a - k);
		da = ba + aa * d0;
		if (da == 0.0) {
			da = e2;
		}
		ba = (1.0 + k * 2.0) - a + z;
		aa = k * (a - k);
		ca = ba + aa / c0;
		if (ca == 0.0) {
			ca = e2;
		}
		da = 1.0 / da;
		de = ca * da;
		fa = fa * de;
		d0 = da;
		c0 = ca;
		k  = k + 1.0;
	} while (mc_fabs(de - 1.0) > e1);
	return fa;
}

MC_TARGET_PROC long double mc_igamma_cfracl_approx1(const long double a, const long double z)
{
	const long double e1 = MCLIMITS_EPSILONL;
	const long double e2 = e1 * 2.0L;
	long double fa       = 1.0L - a + z;
	if (fa == 0.0L) {
		fa = e2;
	}
	long double aa, ba, ca, c0 = fa, da, d0 = 0.0L, de = 0.0L;
	long double k = 1.0L;
	do {
		ba = (1.0L + k * 2.0L) - a + z;
		aa = k * (a - k);
		da = ba + aa * d0;
		if (da == 0.0L) {
			da = e2;
		}
		ba = (1.0L + k * 2.0L) - a + z;
		aa = k * (a - k);
		ca = ba + aa / c0;
		if (ca == 0.0L) {
			ca = e2;
		}
		da = 1.0L / da;
		de = ca * da;
		fa = fa * de;
		d0 = da;
		c0 = ca;
		k  = k + 1.0L;
	} while (mc_fabsl(de - 1.0L) > e1);
	return fa;
}

#pragma mark - mc_igamma_p_approx1 -

MC_TARGET_PROC float mc_igamma_pf_approx1(const float a, const float z)
{
//!# Normalised lower incomplete gamma function of a and z.
	float p = MCK_NAN;
	if (!(z <= 0.0f || a < 0.0f)) {
		const float y = a * mc_logf(z) - z;
		if (y < -FLT_MAX_10_EXP) {
			return 0.0f;
		}
		const float w = y - mc_gammalnf_approx0(a, MC_NULLPTR);
		if (z < (a + 1.0f)) {
			p = mc_igamma_pseriesf_approx1(a, z);
			p = p * mc_expf(w);
		} else {
			p = 1.0f / mc_igamma_cfracf_approx1(a, z);
			p = 1.0f - (mc_expf(w) * p);
		}
	}
	return p;
}

MC_TARGET_PROC double mc_igamma_p_approx1(const double a, const double z)
{
//!# Normalised lower incomplete gamma function of a and z.
	double p = MCK_NAN;
	if (!(z <= 0.0 || a < 0.0)) {
		const double y = a * mc_log(z) - z;
		if (y < -DBL_MAX_10_EXP) {
			return 0.0;
		}
		const double w = y - mc_gammaln_approx0(a, MC_NULLPTR);
		if (z < (a + 1.0)) {
			p = mc_igamma_pseries_approx1(a, z);
			p = p * mc_exp(w);
		} else {
			p = 1.0 / mc_igamma_cfrac_approx1(a, z);
			p = 1.0 - (mc_exp(w) * p);
		}
	}
	return p;
}

MC_TARGET_PROC long double mc_igamma_pl_approx1(const long double a, const long double z)
{
//!# Normalised lower incomplete gamma function of a and z.
	long double p = MCK_NAN;
	if (!(z <= 0.0L || a < 0.0L)) {
		const long double y = a * mc_logl(z) - z;
		if (y < -LDBL_MAX_10_EXP) {
			return 0.0L;
		}
		const long double w = y - mc_gammalnl_approx0(a, MC_NULLPTR);
		if (z < (a + 1.0L)) {
			p = mc_igamma_pseriesl_approx1(a, z);
			p = p * mc_expl(w);
		} else {
			p = 1.0L / mc_igamma_cfracl_approx1(a, z);
			p = 1.0L - (mc_expl(w) * p);
		}
	}
	return p;
}

#pragma mark - mc_igamma_q_approx1 -

MC_TARGET_PROC float mc_igamma_qf_approx1(const float a, const float z)
{
//!# Normalised upper incomplete gamma function of a and z.
	float q = MCK_NAN;
	if (!(z <= 0.0f || a < 0.0f)) {
		const float y = a * mc_logf(z) - z;
		if (y < -FLT_MAX_10_EXP) {
			return 0.0f;
		}
		const float w = y - mc_gammalnf_approx0(a, MC_NULLPTR);
		if (z < (a + 1.0f)) {
			q = mc_igamma_pseriesf_approx1(a, z);
			q = 1.0f - (q * mc_expf(w));
		} else {
			q = 1.0f / mc_igamma_cfracf_approx1(a, z);
			q = mc_expf(w) * q;
		}
	}
	return q;
}

MC_TARGET_PROC double mc_igamma_q_approx1(const double a, const double z)
{
//!# Normalised upper incomplete gamma function of a and z.
	double q = MCK_NAN;
	if (!(z <= 0.0 || a < 0.0)) {
		const double y = a * mc_log(z) - z;
		if (y < -DBL_MAX_10_EXP) {
			return 0.0;
		}
		const double w = y - mc_gammaln_approx0(a, MC_NULLPTR);
		if (z < (a + 1.0)) {
			q = mc_igamma_pseries_approx1(a, z);
			q = 1.0 - (q * mc_exp(w));
		} else {
			q = 1.0 / mc_igamma_cfrac_approx1(a, z);
			q = mc_exp(w) * q;
		}
	}
	return q;
}

MC_TARGET_PROC long double mc_igamma_ql_approx1(const long double a, const long double z)
{
//!# Normalised upper incomplete gamma function of a and z.
	long double q = MCK_NAN;
	if (!(z <= 0.0L || a < 0.0L)) {
		const long double y = a * mc_logl(z) - z;
		if (y < -LDBL_MAX_10_EXP) {
			return 0.0L;
		}
		const long double w = y - mc_gammalnl_approx0(a, MC_NULLPTR);
		if (z < (a + 1.0L)) {
			q = mc_igamma_pseriesl_approx1(a, z);
			q = 1.0L - (q * mc_expl(w));
		} else {
			q = 1.0L / mc_igamma_cfracl_approx1(a, z);
			q = mc_expl(w) * q;
		}
	}
	return q;
}

#pragma mark - mc_igamma_p_taylor_approx2 -

MC_TARGET_PROC float mc_igamma_p_taylorf_approx2(const float a, const float z)
{
	const float x = a * mc_expf(-z) * (1.0f / mc_gammaf(a + 1.0f));
	const float d = mc_powf(z, x);
	float b, c, p;
	if (a <= 0.0f || z <= 0.0f) {
		return 0.0f;
	}
	b = a;
	c = 1.0f;
	p = 1.0f;
	do {
		b = b + 1.0f;
		c = c * z * (1.0f / b);
		p = p + c;
	} while (c < p * MCLIMITS_EPSILONF);
	return (p * d) * (1.0f / a);
}

MC_TARGET_PROC double mc_igamma_p_taylor_approx2(const double a, const double z)
{
	const double x = a * mc_exp(-z) * (1.0 / mc_gamma(a + 1.0));
	const double d = mc_pow(z, x);
	double b, c, p;
	if (a <= 0.0 || z <= 0.0) {
		return 0.0;
	}
	b = a;
	c = 1.0;
	p = 1.0;
	do {
		b = b + 1.0;
		c = c * z * (1.0 / b);
		p = p + c;
	} while (c < p * MCLIMITS_EPSILON);
	return (p * d) * (1.0 / a);
}

MC_TARGET_PROC long double mc_igamma_p_taylorl_approx2(const long double a, const long double z)
{
	const long double x = a * mc_expl(-z) * (1.0L / mc_gammal(a + 1.0L));
	const long double d = mc_powl(z, x);
	long double b, c, p;
	if (a <= 0.0L || z <= 0.0L) {
		return 0.0L;
	}
	b = a;
	c = 1.0L;
	p = 1.0L;
	do {
		b = b + 1.0L;
		c = c * z * (1.0L / b);
		p = p + c;
	} while (c < p * MCLIMITS_EPSILONL);
	return (p * d) * (1.0L / a);
}

#pragma mark - mc_igamma_q_taylor_approx2 -

MC_TARGET_PROC float mc_igamma_q_taylorf_approx2(const float a, const float z)
{
	return 1.0f - mc_igamma_p_taylorf_approx2(a, z);
}

MC_TARGET_PROC double mc_igamma_q_taylor_approx2(const double a, const double z)
{
	return 1.0 - mc_igamma_p_taylor_approx2(a, z);
}

MC_TARGET_PROC long double mc_igamma_q_taylorl_approx2(const long double a, const long  double z)
{
	return 1.0L - mc_igamma_p_taylorl_approx2(a, z);
}

#pragma mark - mc_igamma_q_intarg_approx2 -

MC_TARGET_PROC float mc_igamma_q_intargf_approx2 (const float a, const float z)
{
//!# Normalised upper incomplete gamma function of a and z for a integral.
	float q = mc_expf(-z), w, n;
	if (q > 0.0f && a >= 2.0f) {
		w = q;
		n = 1.0f;
		do {
			w = (w / n) * z;
			q = q + w;
			n = n + 1.0f;
		} while (n < a);
	}
	return q;
}

MC_TARGET_PROC double mc_igamma_q_intarg_approx2(const double a, const double z)
{
//!# Normalised upper incomplete gamma function of a and z for a integral.
	double q = mc_exp(-z), w, n;
	if (q > 0.0 && a >= 2.0) {
		w = q;
		n = 1.0;
		do {
			w = (w / n) * z;
			q = q + w;
			n = n + 1.0;
		} while (n < a);
	}
	return q;
}

MC_TARGET_PROC long double mc_igamma_q_intargl_approx2(const long double a, const long double z)
{
//!# Normalised upper incomplete gamma function of a and z for a integral.
	long double q = mc_expl(-z), w, n;
	if (q > 0.0L && a >= 2.0L) {
		w = q;
		n = 1.0L;
		do {
			w = (w / n) * z;
			q = q + w;
			n = n + 1.0L;
		} while (n < a);
	}
	return q;
}

#pragma mark - mc_igamma_q_halff_approx2 -

MC_TARGET_PROC float mc_igamma_q_halff_approx2(const float a, const float z)
{
//!# Normalised upper incomplete gamma function of a and z for a half-integral.
	float q = mc_erfcf(mc_sqrtf(z)), s, w, n;
	if (q != 0.0f && a > 1.0f) {
		w = 2.0f * mc_expf(-z) * mc_sqrtf(z / MCK_KF(MCK_PI));
		s = w;
		n = 2.0f;
		do {
			w = (w / (n - 0.5f)) * z;
			s = s + w;
			n = n + 1.0f;
		} while (n < a);
		q = q + s;
	}
	return q;
}

MC_TARGET_PROC double mc_igamma_q_half_approx2(const double a, const double z)
{
//!# Normalised upper incomplete gamma function of a and z for a half-integral.
	double q = mc_erfc(mc_sqrt(z)), s, w, n;
	if (q != 0.0 && a > 1.0) {
		w = 2.0 * mc_exp(-z) * mc_sqrt(z / MCK_K(MCK_PI));
		s = w;
		n = 2.0;
		do {
			w = (w / (n - 0.5)) * z;
			s = s + w;
			n = n + 1.0;
		} while (n < a);
		q = q + s;
	}
	return q;
}

MC_TARGET_PROC long double mc_igamma_q_halfl_approx2(const long double a, const long double z)
{
//!# Normalised upper incomplete gamma function of a and z for a half-integral.
	long double q = mc_erfcl(mc_sqrtl(z)), s, w, n;
	if (q != 0.0L && a > 1.0L) {
		w = 2.0L * mc_expl(-z) * mc_sqrtl(z / MCK_KL(MCK_PI));
		s = w;
		n = 2.0L;
		do {
			w = (w / (n - 0.5L)) * z;
			s = s + w;
			n = n + 1.0L;
		} while (n < a);
		q = q + s;
	}
	return q;
}

#pragma mark - mc_igamma_p_approx2 -

MC_TARGET_PROC float mc_igamma_pf_approx2(const float a, const float z)
{
	float p = MCK_INFP;
	if (z <= 0.0f || a < 0.0f) {
		p = MCK_NAN;
	} else if (a == 0.0f && z == 0.0f) {
		p = 0.5f;
	} else if (z == 0.0f) {
		p = 0.0f;
	} else if (a == 0.0f) {
		p = 1.0f;
	} else {
		if (mc_fisintf(2.0f * a)) {
			if (a == 0.5f) {
				p = mc_sqrtf(z);
				return p < 1.0f ? mc_erff(p) : 1.0f - mc_erfcf(p);
			} else if (a == 1.0f) {
				return -mc_expm1f(-z);
			} else if (a < 30.0f && a < (z + 1.0f)) {
				if (mc_fisintf(a) && z > 0.6f) {
					return 1.0f - mc_igamma_q_intargf_approx2(a, z);
				} else if (z > 0.2f) {
					return 1.0f - mc_igamma_q_halff_approx2(a, z);
				}
			}
		}
		p = mc_igamma_pf_approx1(a, z);
	}
	return p;
}

MC_TARGET_PROC double mc_igamma_p_approx2(const double a, const double z)
{
	double p = MCK_INFP;
	if (z <= 0.0 || a < 0.0) {
		p = MCK_NAN;
	} else if (a == 0.0 && z == 0.0) {
		p = 0.5;
	} else if (z == 0.0) {
		p = 0.0;
	} else if (a == 0.0) {
		p = 1.0;
	} else {
		if (mc_fisint(2.0 * a)) {
			if (a == 0.5) {
				p = mc_sqrt(z);
				return p < 1.0 ? mc_erf(p) : 1.0 - mc_erfc(p);
			} else if (a == 1.0) {
				return -mc_expm1(-z);
			} else if (a < 30.0 && a < (z + 1.0)) {
				if (mc_fisint(a) && z > 0.6) {
					return 1.0 - mc_igamma_q_intarg_approx2(a, z);
				} else if (z > 0.2) {
					return 1.0 - mc_igamma_q_half_approx2(a, z);
				}
			}
		}
		p = mc_igamma_p_approx1(a, z);
	}
	return p;
}

MC_TARGET_PROC long double mc_igamma_pl_approx2(const long double a, const long double z)
{
	long double p = MCK_INFP;
	if (z <= 0.0L || a < 0.0L) {
		p = MCK_NAN;
	} else if (a == 0.0L && z == 0.0L) {
		p = 0.5L;
	} else if (z == 0.0L) {
		p = 0.0L;
	} else if (a == 0.0L) {
		p = 1.0L;
	} else {
		if (mc_fisintl(2.0L * a)) {
			if (a == 0.5L) {
				p = mc_sqrtl(z);
				return p < 1.0L ? mc_erfl(p) : 1.0L - mc_erfcl(p);
			} else if (a == 1.0L) {
				return -mc_expm1l(-z);
			} else if (a < 30.0L && a < (z + 1.0L)) {
				if (mc_fisintl(a) && z > 0.6L) {
					return 1.0L - mc_igamma_q_intargl_approx2(a, z);
				} else if (z > 0.2L) {
					return 1.0L - mc_igamma_q_halfl_approx2(a, z);
				}
			}
		}
		p = mc_igamma_pl_approx1(a, z);
	}
	return p;
}

#pragma mark - mc_gamma_q_approx2 -

MC_TARGET_PROC float mc_igamma_qf_approx2(const float a, const float z)
{
	float q = MCK_INFP;
	if (z <= 0.0f || a < 0.0f) {
		q = MCK_NAN;;
	} else if (a == 0.0f && z == 0.0f) {
		q = 0.5f;
	} else if (z == 0.0f) {
		q = 1.0f;
	} else if (a == 0.0f) {
		q = 0.0f;
	} else {
		if (mc_fisintf(2.0f * a)) {
			if (a == 0.5f) {
				q = mc_sqrtf(z);
				return q < 1.0f ? 1.0f - mc_erff(q) : mc_erfcf(q);
			} else if (a == 1.0f) {
				return mc_expf(-z);
			} else if (a < 30.0f && a < (z + 1.0f)) {
				if (mc_fisintf(a) && z > 0.6f) {
					return mc_igamma_q_intargf_approx2(a, z);
				} else if (z > 0.2f) {
					return mc_igamma_q_halff_approx2(a, z);
				}
			}
		}
		q = mc_igamma_qf_approx1(a, z);
	}
	return q;
}

MC_TARGET_PROC double mc_igamma_q_approx2(const double a, const double z)
{
	double q = MCK_INFP;
	if (z <= 0.0 || a < 0.0) {
		q = MCK_NAN;
	} else if (a == 0.0 && z == 0.0) {
		q = 0.5;
	} else if (z == 0.0) {
		q = 1.0;
	} else if (a == 0.0) {
		q = 0.0;
	} else {
		if (mc_fisint(2.0 * a)) {
			if (a == 0.5) {
				q = mc_sqrt(z);
				return q < 1.0 ? 1.0 - mc_erf(q) : mc_erfc(q);
			} else if (a == 1.0) {
				return mc_exp(-z);
			} else if (a < 30.0 && a < (z + 1.0)) {
				if (mc_fisint(a) && z > 0.6) {
					return mc_igamma_q_intarg_approx2(a, z);
				} else if (z > 0.2) {
					return mc_igamma_q_half_approx2(a, z);
				}
			}
		}
		q = mc_igamma_q_approx1(a, z);
	}
	return q;
}

MC_TARGET_PROC long double mc_igamma_ql_approx2(const long double a, const long double z)
{
	long double q = MCK_INFP;
	if (z <= 0.0L || a < 0.0L) {
		q = MCK_NAN;
	} else if (a == 0.0L && z == 0.0L) {
		q = 0.5L;
	} else if (z == 0.0L) {
		q = 1.0L;
	} else if (a == 0.0L) {
		q = 0.0L;
	} else {
		if (mc_fisintl(2.0L * a)) {
			if (a == 0.5L) {
				q = mc_sqrtl(z);
				return q < 1.0L ? 1.0L - mc_erfl(q) : mc_erfcl(q);
			} else if (a == 1.0L) {
				return mc_expl(-z);
			} else if (a < 30.0L && a < (z + 1.0L)) {
				if (mc_fisintl(a) && z > 0.6L) {
					return mc_igamma_q_intargl_approx2(a, z);
				} else if (z > 0.2L) {
					return mc_igamma_q_halfl_approx2(a, z);
				}
			}
		}
		q = mc_igamma_ql_approx1(a, z);
	}
	return q;
}

#endif /* !MC_IGAMMA_H */

/* EOF */