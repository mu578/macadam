//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logcf.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>

#ifndef MC_LOGCF_H
#define MC_LOGCF_H

#pragma mark - mc_logcf -

MC_TARGET_FUNC float mc_logcff(const float x, const float i, const float d, const float tol)
{
//!# Computing a continued fraction approximation to the series:
//!# 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... + x^n/(i+n*d).
//!# Scale factor: huge = 2^64, tiny = 1/2^64.
	const float huge = +1.84467440737095516160000000000000000000E+19f;
	const float tiny = +5.42101086242752217003726400434970855713E-20f;

	float r = 0.0f, a1, a2, b1, b2, c1, c2, c3, c4;

	if (x != 0.0f && d >= 0.0f && d >= 0.0f) {
		c1 = 2.0f * d;
		c2 = i + d;
		c4 = c2 + d;
		a1 = c2;
		b1 = i*(c2 - i*x);
		b2 = d * d * x;
		a2 = c4 * c2 - b2;
		b2 = c4 * b1 - i * b2;

		do {
			c3 = c2 * c2 * x;
			c2 = c2 + d;
			c4 = c4 + d;
			a1 = c4 * a2 - c3 * a1;
			b1 = c4 * b2 - c3 * b1;
			c3 = c1 * c1 * x;
			c1 = c1 + d;
			c4 = c4 + d;
			a2 = c4 * a1 - c3 * a2;
			b2 = c4 * b1 - c3 * b2;
			c3 = mc_fabsf(b2);
			if (c3 > huge) {
				a1 = a1 * tiny;
				b1 = b1 * tiny;
				a2 = a2 * tiny;
				b2 = b2 * tiny;
			} else if (c3 < tiny) {
				a1 = a1 * huge;
				b1 = b1 * huge;
				a2 = a2 * huge;
				b2 = b2 * huge;
			}
		} while (mc_fabsf(a2 * b1 - a1 * b2) > mc_fabsf(tol * b1 * b2));
		r = a2 / b2;
	}
	return r;
}

MC_TARGET_FUNC double mc_logcf(const double x, const double i, const double d, const double tol)
{
//!# Computing a continued fraction approximation to the series:
//!# 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... + x^n/(i+n*d).
//!# Scale factor: huge = 2^256, tiny = 1/2^256.
	const double huge = +1.1579208923731619542357098500868790785327E+77;
	const double tiny = +8.6361685550944446253863518628003995711160E-78;

	double r = 0.0, a1, a2, b1, b2, c1, c2, c3, c4;

	if (x != 0.0 && d >= 0.0 && d >= 0.0) {
		c1 = 2.0 * d;
		c2 = i + d;
		c4 = c2 + d;
		a1 = c2;
		b1 = i*(c2 - i*x);
		b2 = d * d * x;
		a2 = c4 * c2 - b2;
		b2 = c4 * b1 - i * b2;

		do {
			c3 = c2 * c2 * x;
			c2 = c2 + d;
			c4 = c4 + d;
			a1 = c4 * a2 - c3 * a1;
			b1 = c4 * b2 - c3 * b1;
			c3 = c1 * c1 * x;
			c1 = c1 + d;
			c4 = c4 + d;
			a2 = c4 * a1 - c3 * a2;
			b2 = c4 * b1 - c3 * b2;
			c3 = mc_fabs(b2);
			if (c3 > huge) {
				a1 = a1 * tiny;
				b1 = b1 * tiny;
				a2 = a2 * tiny;
				b2 = b2 * tiny;
			} else if (c3 < tiny) {
				a1 = a1 * huge;
				b1 = b1 * huge;
				a2 = a2 * huge;
				b2 = b2 * huge;
			}
		} while (mc_fabs(a2 * b1 - a1 * b2) > mc_fabs(tol * b1 * b2));
		r = a2 / b2;
	}
	return r;
}

MC_TARGET_FUNC long double mc_logcfl(const long double x, const long double i, const long double d, const long double tol)
{
//!# Computing a continued fraction approximation to the series:
//!# 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... + x^n/(i+n*d).
//!# Scale factor: huge = 2^1024, tiny = 1/2^1024.
#	if MC_TARGET_HAVE_LONG_DOUBLE
	const long double huge = +1.797693134862315907729305190789024733617976978942306572734301E+308L;
	const long double tiny = +5.562684646268003457725581793331010160548039951155829576383319E-309L;
#	else
	const long double huge = +1.1579208923731619542357098500868790785327E+77L;
	const long double tiny = +8.6361685550944446253863518628003995711160E-78L;
#	endif

	long double r = 0.0L, a1, a2, b1, b2, c1, c2, c3, c4;

	if (x != 0.0L && d >= 0.0L && d >= 0.0L) {
		c1 = 2.0L * d;
		c2 = i + d;
		c4 = c2 + d;
		a1 = c2;
		b1 = i*(c2 - i*x);
		b2 = d * d * x;
		a2 = c4 * c2 - b2;
		b2 = c4 * b1 - i * b2;

		do {
			c3 = c2 * c2 * x;
			c2 = c2 + d;
			c4 = c4 + d;
			a1 = c4 * a2 - c3 * a1;
			b1 = c4 * b2 - c3 * b1;
			c3 = c1 * c1 * x;
			c1 = c1 + d;
			c4 = c4 + d;
			a2 = c4 * a1 - c3 * a2;
			b2 = c4 * b1 - c3 * b2;
			c3 = mc_fabsl(b2);
			if (c3 > huge) {
				a1 = a1 * tiny;
				b1 = b1 * tiny;
				a2 = a2 * tiny;
				b2 = b2 * tiny;
			} else if (c3 < tiny) {
				a1 = a1 * huge;
				b1 = b1 * huge;
				a2 = a2 * huge;
				b2 = b2 * huge;
			}
		} while (mc_fabsl(a2 * b1 - a1 * b2) > mc_fabsl(tol * b1 * b2));
		r = a2 / b2;
	}
	return r;
}

#endif /* !MC_LOGCF_H */

/* EOF */