//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lartgp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_LAPACKE_LARTGP_H
#define MC_LAPACKE_LARTGP_H

#pragma mark - mc_lapack_slartgp -

MC_TARGET_FUNC void mc_lapack_slartgp(float f, float g, float * cs, float * sn, float * r)
{
	const float one = 1.0f, zero = 0.0f;

	const float safmn2 = mc_lapack_slamch('W');
	const float safmx2 = one / safmn2;

	int count, i;
	float f1, g1, scale;

	if (g == zero) {
		*cs = mc_copysignf(one, f);
		*sn = zero;
		*r  = mc_fabsf(f);
	} else if (f == zero) {
		*cs = zero;
		*sn = mc_copysignf(one, g);
		*r  = mc_fabsf(g);
	} else {
		f1    = f;
		g1    = g;
		scale = mc_fmaxf(mc_fabsf(f1), mc_fabsf(g1));
		if (scale >= safmx2) {
			count = 0;
F10:
			count = count + 1;
			f1    = f1 * safmn2;
			g1    = g1 * safmn2;
			scale = mc_fmaxf(mc_fabsf(f1), mc_fabsf(g1));
			if (scale >= safmx2) {
				goto F10;
			}
			*r  = mc_sqrtf(mc_raise2f(f1) + mc_raise2f(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
			for (i = 1; i <= count; ++i) {
				*r = (*r) * safmx2;
			}
		} else if (scale <= safmn2) {
			count = 0;
F30:
			count = count + 1;
			f1    = f1 * safmx2;
			g1    = g1 * safmx2;
			scale = mc_fmaxf(mc_fabsf(f1), mc_fabsf(g1));
			if (scale <= safmn2) {
				goto F30;
			}
			*r  = mc_sqrtf(mc_raise2f(f1) + mc_raise2f(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
			for (i = 1; i <= count; ++i) {
				*r = (*r) * safmn2;
			}
		} else {
			*r  = mc_sqrtf(mc_raise2f(f1) + mc_raise2f(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
		}
		if ((*r) < zero) {
			*cs = -(*cs);
			*sn = -(*sn);
			*r  = -(*r);
		}
	}
}

#pragma mark - mc_lapack_dlartgp -

MC_TARGET_FUNC void mc_lapack_dlartgp(double f, double g, double * cs, double * sn, double * r)
{
	const double one = 1.0, zero = 0.0;

	const double safmn2 = mc_lapack_dlamch('W');
	const double safmx2 = one / safmn2;

	int count, i;
	double f1, g1, scale;

	if (g == zero) {
		*cs = mc_copysign(one, f);
		*sn = zero;
		*r  = mc_fabs(f);
	} else if (f == zero) {
		*cs = zero;
		*sn = mc_copysign(one, g);
		*r  = mc_fabs(g);
	} else {
		f1    = f;
		g1    = g;
		scale = mc_fmax(mc_fabs(f1), mc_fabs(g1));
		if (scale >= safmx2) {
			count = 0;
F10:
			count = count + 1;
			f1    = f1 * safmn2;
			g1    = g1 * safmn2;
			scale = mc_fmax(mc_fabs(f1), mc_fabs(g1));
			if (scale >= safmx2) {
				goto F10;
			}
			*r  = mc_sqrt(mc_raise2(f1) + mc_raise2(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
			for (i = 1; i <= count; ++i) {
				*r = (*r) * safmx2;
			}
		} else if (scale <= safmn2) {
			count = 0;
F30:
			count = count + 1;
			f1    = f1 * safmx2;
			g1    = g1 * safmx2;
			scale = mc_fmax(mc_fabs(f1), mc_fabs(g1));
			if (scale <= safmn2) {
				goto F30;
			}
			*r  = mc_sqrt(mc_raise2(f1) + mc_raise2(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
			for (i = 1; i <= count; ++i) {
				*r = (*r) * safmn2;
			}
		} else {
			*r  = mc_sqrt(mc_raise2(f1) + mc_raise2(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
		}
		if ((*r) < zero) {
			*cs = -(*cs);
			*sn = -(*sn);
			*r  = -(*r);
		}
	}
}

#pragma mark - mc_lapack_llartgp -

MC_TARGET_FUNC void mc_lapack_llartgp(long double f, long double g, long double * cs, long double * sn, long double * r)
{
	const long double one = 1.0L, zero = 0.0L;

	const long double safmn2 = mc_lapack_llamch('W');
	const long double safmx2 = one / safmn2;

	int count, i;
	long double f1, g1, scale;

	if (g == zero) {
		*cs = mc_copysignl(one, f);
		*sn = zero;
		*r  = mc_fabsl(f);
	} else if (f == zero) {
		*cs = zero;
		*sn = mc_copysignl(one, g);
		*r  = mc_fabsl(g);
	} else {
		f1    = f;
		g1    = g;
		scale = mc_fmaxl(mc_fabsl(f1), mc_fabsl(g1));
		if (scale >= safmx2) {
			count = 0;
F10:
			count = count + 1;
			f1    = f1 * safmn2;
			g1    = g1 * safmn2;
			scale = mc_fmaxl(mc_fabsl(f1), mc_fabsl(g1));
			if (scale >= safmx2) {
				goto F10;
			}
			*r  = mc_sqrtl(mc_raise2l(f1) + mc_raise2l(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
			for (i = 1; i <= count; ++i) {
				*r = (*r) * safmx2;
			}
		} else if (scale <= safmn2) {
			count = 0;
F30:
			count = count + 1;
			f1    = f1 * safmx2;
			g1    = g1 * safmx2;
			scale = mc_fmaxl(mc_fabsl(f1), mc_fabsl(g1));
			if (scale <= safmn2) {
				goto F30;
			}
			*r  = mc_sqrtl(mc_raise2l(f1) + mc_raise2l(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
			for (i = 1; i <= count; ++i) {
				*r = (*r) * safmn2;
			}
		} else {
			*r  = mc_sqrtl(mc_raise2l(f1) + mc_raise2l(g1));
			*cs = f1 / (*r);
			*sn = g1 / (*r);
		}
		if ((*r) < zero) {
			*cs = -(*cs);
			*sn = -(*sn);
			*r  = -(*r);
		}
	}
}

#endif /* !MC_LAPACKE_LARTGP_H */

/* EOF */