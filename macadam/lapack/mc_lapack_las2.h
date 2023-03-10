//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_las2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#	include <macadam/lapack/mc_blas.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_fmin.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_LAPACKE_LAS2_H
#define MC_LAPACKE_LAS2_H

#pragma mark - mc_lapack_slas2 -

MC_TARGET_FUNC void mc_lapack_slas2(float f, float g, float h, float * ssmin, float * ssmax)
{
	const float two = 2.0f, one = 1.0f, zero = 0.0f;

	float c, fa, ga, ha, as, at, au, fhmn, fhmx;

	fa   = mc_fabsf(f);
	ga   = mc_fabsf(g);
	ha   = mc_fabsf(h);
	fhmn = mc_fminf(fa, ha);
	fhmx = mc_fmaxf(fa, ha);
	if (fhmn == zero) {
		*ssmin = zero;
		if (fhmx == zero) {
			*ssmax = ga;
		} else {
			*ssmax = mc_fmaxf(fhmx,ga) * mc_sqrtf(one + mc_raise2f(mc_fminf(fhmx, ga) / mc_fmaxf(fhmx, ga)));
		}
	} else {
		if (ga < fhmx) {
			 as    = fhmn / fhmx + one;
			 at    = (fhmx - fhmn) / fhmx;
			 au    = mc_raise2f(ga / fhmx);
			 c     = two / (mc_sqrtf(mc_raise2f(as) + au) + mc_sqrtf(mc_raise2f(at) + au));
			*ssmin = fhmn * c;
			*ssmax = fhmx / c;
		} else {
			au = fhmx / ga;
			if (au == zero) {
				*ssmin = fhmn * fhmx / ga;
				*ssmax = ga;
			} else {
				 as    = fhmn / fhmx + one;
				 at    = (fhmx - fhmn) / fhmx;
				 c     = one / (mc_sqrtf(one + mc_raise2f(as * au)) + mc_sqrtf(one + mc_raise2f(at * au)));
				*ssmin = fhmn * c * au;
				*ssmin = (*ssmin) + (*ssmin);
				*ssmax = ga / (c + c);
			}
		}
	}
}

#pragma mark - mc_lapack_dlas2 -

MC_TARGET_FUNC void mc_lapack_dlas2(double f, double g, double h, double * ssmin, double * ssmax)
{
	const double two = 2.0, one = 1.0, zero = 0.0;

	double c, fa, ga, ha, as, at, au, fhmn, fhmx;

	fa   = mc_fabs(f);
	ga   = mc_fabs(g);
	ha   = mc_fabs(h);
	fhmn = mc_fmin(fa, ha);
	fhmx = mc_fmax(fa, ha);
	if (fhmn == zero) {
		*ssmin = zero;
		if (fhmx == zero) {
			*ssmax = ga;
		} else {
			*ssmax = mc_fmax(fhmx,ga) * mc_sqrt(one + mc_raise2(mc_fmin(fhmx, ga) / mc_fmax(fhmx, ga)));
		}
	} else {
		if (ga < fhmx) {
			 as    = fhmn / fhmx + one;
			 at    = (fhmx - fhmn) / fhmx;
			 au    = mc_raise2(ga / fhmx);
			 c     = two / (mc_sqrt(mc_raise2(as) + au) + mc_sqrt(mc_raise2(at) + au));
			*ssmin = fhmn * c;
			*ssmax = fhmx / c;
		} else {
			au = fhmx / ga;
			if (au == zero) {
				*ssmin = fhmn * fhmx / ga;
				*ssmax = ga;
			} else {
				 as    = fhmn / fhmx + one;
				 at    = (fhmx - fhmn) / fhmx;
				 c     = one / (mc_sqrt(one + mc_raise2(as * au)) + mc_sqrt(one + mc_raise2(at * au)));
				*ssmin = fhmn * c * au;
				*ssmin = (*ssmin) + (*ssmin);
				*ssmax = ga / (c + c);
			}
		}
	}
}

#pragma mark - mc_lapack_llas2 -

MC_TARGET_FUNC void mc_lapack_llas2(long double f, long double g, long double h, long double * ssmin, long double * ssmax)
{
	const long double two = 2.0L, one = 1.0L, zero = 0.0L;

	long double c, fa, ga, ha, as, at, au, fhmn, fhmx;

	fa   = mc_fabsl(f);
	ga   = mc_fabsl(g);
	ha   = mc_fabsl(h);
	fhmn = mc_fminl(fa, ha);
	fhmx = mc_fmaxl(fa, ha);
	if (fhmn == zero) {
		*ssmin = zero;
		if (fhmx == zero) {
			*ssmax = ga;
		} else {
			*ssmax = mc_fmaxl(fhmx,ga) * mc_sqrtl(one + mc_raise2l(mc_fminl(fhmx, ga) / mc_fmaxl(fhmx, ga)));
		}
	} else {
		if (ga < fhmx) {
			 as    = fhmn / fhmx + one;
			 at    = (fhmx - fhmn) / fhmx;
			 au    = mc_raise2l(ga / fhmx);
			 c     = two / (mc_sqrtl(mc_raise2l(as) + au) + mc_sqrtl(mc_raise2l(at) + au));
			*ssmin = fhmn * c;
			*ssmax = fhmx / c;
		} else {
			au = fhmx / ga;
			if (au == zero) {
				*ssmin = fhmn * fhmx / ga;
				*ssmax = ga;
			} else {
				 as    = fhmn / fhmx + one;
				 at    = (fhmx - fhmn) / fhmx;
				 c     = one / (mc_sqrtl(one + mc_raise2l(as * au)) + mc_sqrtl(one + mc_raise2l(at * au)));
				*ssmin = fhmn * c * au;
				*ssmin = (*ssmin) + (*ssmin);
				*ssmax = ga / (c + c);
			}
		}
	}
}

#endif /* !MC_LAPACKE_LAS2_H */

/* EOF */