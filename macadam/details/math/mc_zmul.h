//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zmul.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>

#ifndef MC_ZMUL_H
#define MC_ZMUL_H

#pragma mark - mc_zmul -

MC_TARGET_PROC void mc_zmulf(float * c_r, float * c_i
	, float a_r, float a_i
	, float b_r, float b_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = (a_r * b_r) - (a_i * b_i);
	*c_i = (a_r * b_i) + (a_i * b_r);
#	else
	float a = a_r * b_r, b = a_i * b_i;
	float c = a_r * b_i, d = a_i * b_r;
	int w   = 0;
	*c_r    = a - b;
	*c_i    = c + d;
	if (mc_isnan(*c_r) && mc_isnan(*c_i)) {
		if (mc_isinf(a_r) || mc_isinf(a_i)) {
			a_r = mc_copysignf(mc_isinf(a_r) ? 1.0f : 0.0f, a_r);
			a_i = mc_copysignf(mc_isinf(a_i) ? 1.0f : 0.0f, a_i);
			if (mc_isnan(b_r)) {
				b_r = mc_copysignf(0.0f, b_r);
			}
			if (mc_isnan(b_i)) {
				b_i = mc_copysignf(0.0f, b_i);
			}
			w = 1;
		}
		if (mc_isinf(b_r) || mc_isinf(b_i)) {
			b_r = mc_copysignf(mc_isinf(b_r) ? 1.0f : 0.0f, b_r);
			b_i = mc_copysignf(mc_isinf(b_i) ? 1.0f : 0.0f, b_i);
			if (mc_isnan(a_r)) {
				a_r = mc_copysignf(0.0f, a_r);
			}
			if (mc_isnan(a_i)) {
				a_i = mc_copysignf(0.0f, a_i);
			}
			w = 1;
		}
		if (w == 0 && (mc_isinf(a) || mc_isinf(b) || mc_isinf(c) || mc_isinf(d))) {
			if (mc_isnan(a_r)) {
				a_r = mc_copysignf(0.0f, a_r);
			}
			if (mc_isnan(a_i)) {
				a_i = mc_copysignf(0.0f, a_i);
			}
			if (mc_isnan(b_r)) {
				b_r = mc_copysignf(0.0f, b_r);
			}
			if (mc_isnan(b_i)) {
				b_i = mc_copysignf(0.0f, b_i);
			}
			w = 1;
		}
		if (w) {
			*c_r = MCK_INFP * (a_r * b_r - a_i * b_i);
			*c_i = MCK_INFP * (a_r * b_i + a_i * b_r);
		}
	}
#	endif
}

MC_TARGET_PROC void mc_zmul(double * c_r, double * c_i
	, double a_r, double a_i
	, double b_r, double b_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = (a_r * b_r) - (a_i * b_i);
	*c_i = (a_r * b_i) + (a_i * b_r);
#	else
	double a = a_r * b_r, b = a_i * b_i;
	double c = a_r * b_i, d = a_i * b_r;
	int w    = 0;
	*c_r     = a - b;
	*c_i     = c + d;
	if (mc_isnan(*c_r) && mc_isnan(*c_i)) {
		if (mc_isinf(a_r) || mc_isinf(a_i)) {
			a_r = mc_copysign(mc_isinf(a_r) ? 1.0 : 0.0, a_r);
			a_i = mc_copysign(mc_isinf(a_i) ? 1.0 : 0.0, a_i);
			if (mc_isnan(b_r)) {
				b_r = mc_copysign(0.0, b_r);
			}
			if (mc_isnan(b_i)) {
				b_i = mc_copysign(0.0, b_i);
			}
			w = 1;
		}
		if (mc_isinf(b_r) || mc_isinf(b_i)) {
			b_r = mc_copysign(mc_isinf(b_r) ? 1.0 : 0.0, b_r);
			b_i = mc_copysign(mc_isinf(b_i) ? 1.0 : 0.0, b_i);
			if (mc_isnan(a_r)) {
				a_r = mc_copysign(0.0, a_r);
			}
			if (mc_isnan(a_i)) {
				a_i = mc_copysign(0.0, a_i);
			}
			w = 1;
		}
		if (w == 0 && (mc_isinf(a) || mc_isinf(b) || mc_isinf(c) || mc_isinf(d))) {
			if (mc_isnan(a_r)) {
				a_r = mc_copysign(0.0, a_r);
			}
			if (mc_isnan(a_i)) {
				a_i = mc_copysign(0.0, a_i);
			}
			if (mc_isnan(b_r)) {
				b_r = mc_copysign(0.0, b_r);
			}
			if (mc_isnan(b_i)) {
				b_i = mc_copysign(0.0, b_i);
			}
			w = 1;
		}
		if (w) {
			*c_r = MCK_INFP * (a_r * b_r - a_i * b_i);
			*c_i = MCK_INFP * (a_r * b_i + a_i * b_r);
		}
	}
#	endif
}

MC_TARGET_PROC void mc_zmull(long double * c_r, long double * c_i
	, long double a_r, long double a_i
	, long double b_r, long double b_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = (a_r * b_r) - (a_i * b_i);
	*c_i = (a_r * b_i) + (a_i * b_r);
#	else
	long double a = a_r * b_r, b = a_i * b_i;
	long double c = a_r * b_i, d = a_i * b_r;
	int w         = 0;
	*c_r          = a - b;
	*c_i          = c + d;
	if (mc_isnan(*c_r) && mc_isnan(*c_i)) {
		if (mc_isinf(a_r) || mc_isinf(a_i)) {
			a_r = mc_copysignl(mc_isinf(a_r) ? 1.0L : 0.0L, a_r);
			a_i = mc_copysignl(mc_isinf(a_i) ? 1.0L : 0.0L, a_i);
			if (mc_isnan(b_r)) {
				b_r = mc_copysignl(0.0L, b_r);
			}
			if (mc_isnan(b_i)) {
				b_i = mc_copysignl(0.0L, b_i);
			}
			w = 1;
		}
		if (mc_isinf(b_r) || mc_isinf(b_i)) {
			b_r = mc_copysignl(mc_isinf(b_r) ? 1.0L : 0.0L, b_r);
			b_i = mc_copysignl(mc_isinf(b_i) ? 1.0L : 0.0L, b_i);
			if (mc_isnan(a_r)) {
				a_r = mc_copysignl(0.0L, a_r);
			}
			if (mc_isnan(a_i)) {
				a_i = mc_copysignl(0.0L, a_i);
			}
			w = 1;
		}
		if (w == 0 && (mc_isinf(a) || mc_isinf(b) || mc_isinf(c) || mc_isinf(d))) {
			if (mc_isnan(a_r)) {
				a_r = mc_copysignl(0.0L, a_r);
			}
			if (mc_isnan(a_i)) {
				a_i = mc_copysignl(0.0L, a_i);
			}
			if (mc_isnan(b_r)) {
				b_r = mc_copysignl(0.0L, b_r);
			}
			if (mc_isnan(b_i)) {
				b_i = mc_copysignl(0.0L, b_i);
			}
			w = 1;
		}
		if (w) {
			*c_r = MCK_INFP * (a_r * b_r - a_i * b_i);
			*c_i = MCK_INFP * (a_r * b_i + a_i * b_r);
		}
	}
#	endif
}

#endif /* !MC_ZMUL_H */

/* EOF */