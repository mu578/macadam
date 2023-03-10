//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasv2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_LAPACKE_LASV2_H
#define MC_LAPACKE_LASV2_H

#pragma mark - mc_lapack_slasv2 -

MC_TARGET_FUNC void mc_lapack_slasv2(float f, float g, float h, float * ssmin, float * ssmax, float * snr, float * csr, float * snl, float * csl)
{
	const float four = 4.0f, two = 2.0f, one = 1.0f, half = 0.5f, zero = 0.0f;

	int gasmal, swap;
	int pmax;
	float a, clt, crt, d, fa, ft, ga, gt, ha, ht, l, m,
	      mm, r, s, slt, srt, t, temp, tsign, tt;

	ft = f;
	fa = mc_fabsf(ft);
	ht = h;
	ha = mc_fabsf(h);

	pmax = 1;
	swap = ha > fa;
	if (swap) {
		pmax = 3;
		temp = ft;
		ft   = ht;
		ht   = temp;
		temp = fa;
		fa   = ha;
		ha   = temp;
	}
	gt = g;
	ga = mc_fabsf(gt);
	if (ga == zero) {
		*ssmin = ha;
		*ssmax = fa;
		 clt   = one;
		 crt   = one;
		 slt   = zero;
		 srt   = zero;
	} else {
		gasmal = 1;
		if (ga > fa) {
			pmax = 2;
			if (fa / ga < mc_lapack_slamch('E')) {
				 gasmal = 0;
				*ssmax  = ga;
				if (ha > one) {
					*ssmin = fa / (ga / ha);
				} else {
					*ssmin = fa / ga * ha;
				}
				clt = one;
				slt = ht / gt;
				srt = one;
				crt = ft / gt;
			}
		}
		if (gasmal) {
			d = fa - ha;
			if (d == fa) {
				l = one;
			} else {
				l = d / fa;
			}
			m  = gt / ft;
			t  = two - l;
			mm = m * m;
			tt = t * t;
			s  = mc_sqrtf(tt + mm);
			if (l == zero) {
				r = mc_fabsf(m);
			} else {
				r = mc_sqrtf(l * l + mm);
			}
			 a     = (s + r) * half;
			*ssmin = ha / a;
			*ssmax = fa * a;
			if (mm == zero) {
				if (l == zero) {
					t = mc_copysignf(two, ft) * mc_copysignf(one, gt);
				} else {
					t = gt / mc_copysignf(d, ft) + m / t;
				}
			} else {
				t = (m / (s + t) + m / (r + l)) * (a + one);
			}
			l   = mc_sqrtf(t * t + four);
			crt = two / l;
			srt = t / l;
			clt = (crt + srt * m) / a;
			slt = ht / ft * srt / a;
		}
	}
	if (swap) {
		*csl = srt;
		*snl = crt;
		*csr = slt;
		*snr = clt;
	} else {
		*csl = clt;
		*snl = slt;
		*csr = crt;
		*snr = srt;
	}

	tsign = zero;
	if (pmax == 1) {
		tsign = mc_copysignf(one, *csr) * mc_copysignf(one, *csl) * mc_copysignf(one, f);
	} else if (pmax == 2) {
		tsign = mc_copysignf(one, *snr) * mc_copysignf(one, *csl) * mc_copysignf(one, g);
	} else if (pmax == 3) {
		tsign = mc_copysignf(one, *snr) * mc_copysignf(one, *snl) * mc_copysignf(one, h);
	}
	*ssmax = mc_copysignf(*ssmax, tsign);
	*ssmin = mc_copysignf(*ssmin, tsign * mc_copysignf(one, f) * mc_copysignf(one, h));
}

#pragma mark - mc_lapack_dlasv2 -

MC_TARGET_FUNC void mc_lapack_dlasv2(double f, double g, double h, double * ssmin, double * ssmax, double * snr, double * csr, double * snl, double * csl)
{
	const double four = 4.0, two = 2.0, one = 1.0, half = 0.5, zero = 0.0;

	int gasmal, swap;
	int pmax;
	double a, clt, crt, d, fa, ft, ga, gt, ha, ht, l, m,
	       mm, r, s, slt, srt, t, temp, tsign, tt;

	ft = f;
	fa = mc_fabs(ft);
	ht = h;
	ha = mc_fabs(h);

	pmax = 1;
	swap = ha > fa;
	if (swap) {
		pmax = 3;
		temp = ft;
		ft   = ht;
		ht   = temp;
		temp = fa;
		fa   = ha;
		ha   = temp;
	}
	gt = g;
	ga = mc_fabs(gt);
	if (ga == zero) {
		*ssmin = ha;
		*ssmax = fa;
		 clt   = one;
		 crt   = one;
		 slt   = zero;
		 srt   = zero;
	} else {
		gasmal = 1;
		if (ga > fa) {
			pmax = 2;
			if (fa / ga < mc_lapack_dlamch('E')) {
				 gasmal = 0;
				*ssmax  = ga;
				if (ha > one) {
					*ssmin = fa / (ga / ha);
				} else {
					*ssmin = fa / ga * ha;
				}
				clt = one;
				slt = ht / gt;
				srt = one;
				crt = ft / gt;
			}
		}
		if (gasmal) {
			d = fa - ha;
			if (d == fa) {
				l = one;
			} else {
				l = d / fa;
			}
			m  = gt / ft;
			t  = two - l;
			mm = m * m;
			tt = t * t;
			s  = mc_sqrt(tt + mm);
			if (l == zero) {
				r = mc_fabs(m);
			} else {
				r = mc_sqrt(l * l + mm);
			}
			 a     = (s + r) * half;
			*ssmin = ha / a;
			*ssmax = fa * a;
			if (mm == zero) {
				if (l == zero) {
					t = mc_copysign(two, ft) * mc_copysign(one, gt);
				} else {
					t = gt / mc_copysign(d, ft) + m / t;
				}
			} else {
				t = (m / (s + t) + m / (r + l)) * (a + one);
			}
			l   = mc_sqrt(t * t + four);
			crt = two / l;
			srt = t / l;
			clt = (crt + srt * m) / a;
			slt = ht / ft * srt / a;
		}
	}
	if (swap) {
		*csl = srt;
		*snl = crt;
		*csr = slt;
		*snr = clt;
	} else {
		*csl = clt;
		*snl = slt;
		*csr = crt;
		*snr = srt;
	}

	tsign = zero;
	if (pmax == 1) {
		tsign = mc_copysign(one, *csr) * mc_copysign(one, *csl) * mc_copysign(one, f);
	} else if (pmax == 2) {
		tsign = mc_copysign(one, *snr) * mc_copysign(one, *csl) * mc_copysign(one, g);
	} else if (pmax == 3) {
		tsign = mc_copysign(one, *snr) * mc_copysign(one, *snl) * mc_copysign(one, h);
	}
	*ssmax = mc_copysign(*ssmax, tsign);
	*ssmin = mc_copysign(*ssmin, tsign * mc_copysign(one, f) * mc_copysign(one, h));
}

#pragma mark - mc_lapack_llasv2 -

MC_TARGET_FUNC void mc_lapack_llasv2(long double f, long double g, long double h, long double * ssmin, long double * ssmax, long double * snr, long double * csr, long double * snl, long double * csl)
{
	const long double four = 4.0L, two = 2.0L, one = 1.0L, half = 0.5L, zero = 0.0L;

	int gasmal, swap;
	int pmax;
	long double a, clt, crt, d, fa, ft, ga, gt, ha, ht, l, m,
	            mm, r, s, slt, srt, t, temp, tsign, tt;

	ft = f;
	fa = mc_fabsl(ft);
	ht = h;
	ha = mc_fabsl(h);

	pmax = 1;
	swap = ha > fa;
	if (swap) {
		pmax = 3;
		temp = ft;
		ft   = ht;
		ht   = temp;
		temp = fa;
		fa   = ha;
		ha   = temp;
	}
	gt = g;
	ga = mc_fabsl(gt);
	if (ga == zero) {
		*ssmin = ha;
		*ssmax = fa;
		 clt   = one;
		 crt   = one;
		 slt   = zero;
		 srt   = zero;
	} else {
		gasmal = 1;
		if (ga > fa) {
			pmax = 2;
			if (fa / ga < mc_lapack_llamch('E')) {
				 gasmal = 0;
				*ssmax  = ga;
				if (ha > one) {
					*ssmin = fa / (ga / ha);
				} else {
					*ssmin = fa / ga * ha;
				}
				clt = one;
				slt = ht / gt;
				srt = one;
				crt = ft / gt;
			}
		}
		if (gasmal) {
			d = fa - ha;
			if (d == fa) {
				l = one;
			} else {
				l = d / fa;
			}
			m  = gt / ft;
			t  = two - l;
			mm = m * m;
			tt = t * t;
			s  = mc_sqrtl(tt + mm);
			if (l == zero) {
				r = mc_fabsl(m);
			} else {
				r = mc_sqrtl(l * l + mm);
			}
			 a     = (s + r) * half;
			*ssmin = ha / a;
			*ssmax = fa * a;
			if (mm == zero) {
				if (l == zero) {
					t = mc_copysignl(two, ft) * mc_copysignl(one, gt);
				} else {
					t = gt / mc_copysignl(d, ft) + m / t;
				}
			} else {
				t = (m / (s + t) + m / (r + l)) * (a + one);
			}
			l   = mc_sqrtl(t * t + four);
			crt = two / l;
			srt = t / l;
			clt = (crt + srt * m) / a;
			slt = ht / ft * srt / a;
		}
	}
	if (swap) {
		*csl = srt;
		*snl = crt;
		*csr = slt;
		*snr = clt;
	} else {
		*csl = clt;
		*snl = slt;
		*csr = crt;
		*snr = srt;
	}

	tsign = zero;
	if (pmax == 1) {
		tsign = mc_copysignl(one, *csr) * mc_copysignl(one, *csl) * mc_copysignl(one, f);
	} else if (pmax == 2) {
		tsign = mc_copysignl(one, *snr) * mc_copysignl(one, *csl) * mc_copysignl(one, g);
	} else if (pmax == 3) {
		tsign = mc_copysignl(one, *snr) * mc_copysignl(one, *snl) * mc_copysignl(one, h);
	}
	*ssmax = mc_copysignl(*ssmax, tsign);
	*ssmin = mc_copysignl(*ssmin, tsign * mc_copysignl(one, f) * mc_copysignl(one, h));
}

#endif /* !MC_LAPACKE_LASV2_H */

/* EOF */