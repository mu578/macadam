//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasq3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_fmin.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_LAPACKE_LASQ4_H
#define MC_LAPACKE_LASQ4_H

#pragma mark - mc_lapack_slasq4 -

MC_TARGET_PROC void mc_lapack_slasq4(const int i0, const int n0, const float * z
	, const int pp
	, const int n0in
	, float     dmin
	, float     dmin1
	, float     dmin2
	, float     dn
	, float     dn1
	, float     dn2
	, float *   tau
	, int *     ttype
	, float *   g
) {
	const float hundrd = 100.0f, two = 2.0f, one = 1.0f, half = 0.5f, third = 0.333f, qurtr = 0.25f, zero = 0.0f;
	const float cnst1 = 0.563f, cnst2 = 1.01f, cnst3 = 1.05f;

	int i4, nn, np;
	float a2, b1, b2, gam, gap1, gap2, s;

	if (dmin <= zero) {
		*tau   = -(dmin);
		*ttype = -1;
		return;
	}

	s  = zero;
	nn = (4 * n0) + pp;
	if (n0in == n0) {
		if (dmin == dn || dmin == dn1) {
			b1 = mc_sqrtf(mc_blas_vector_at(z, nn - 3)) * mc_sqrtf(mc_blas_vector_at(z, nn - 5));
			b2 = mc_sqrtf(mc_blas_vector_at(z, nn - 7)) * mc_sqrtf(mc_blas_vector_at(z, nn - 9));
			a2 = mc_blas_vector_at(z, nn - 7) + mc_blas_vector_at(z, nn - 5);
			if (dmin == dn && dmin1 == dn1) {
				gap2 = dmin2 - a2 - dmin2 * qurtr;
				if (gap2 > zero && gap2 > b2) {
					gap1 = a2 - dn - b2 / gap2 * b2;
				} else {
					gap1 = a2 - dn - (b1 + b2);
				}
				if (gap1 > zero && gap1 > b1) {
					 s     = mc_fmaxf(dn - b1 / gap1 * b1, dmin * half);
					*ttype = -2;
				} else {
					s = zero;
					if (dn > b1) {
						s = dn - b1;
					}
					if (a2 > b1 + b2) {
						s = mc_fminf(s, a2 - (b1 + b2));
					}
					 s     = mc_fmaxf(s, dmin * third);
					*ttype = -3;
				}
			} else {
				*ttype = -4;
				s = dmin * qurtr;
				if (dmin == dn) {
					gam = dn;
					a2  = zero;
					if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
						return;
					}
					b2 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
					np = nn - 9;
				} else {
					np  = nn - (2 * pp);
					gam = dn1;
					if (mc_blas_vector_at(z, np - 4) > mc_blas_vector_at(z, np - 2)) {
						return;
					}
					a2 = mc_blas_vector_at(z, np - 4) / mc_blas_vector_at(z, np - 2);
					if (mc_blas_vector_at(z, nn - 9) > mc_blas_vector_at(z, nn - 11)) {
						return;
					}
					b2 = mc_blas_vector_at(z, nn - 9) / mc_blas_vector_at(z, nn - 11);
					np = nn - 13;
				}
				a2 = a2 + b2;
				for (i4 = np; i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
					if (b2 == zero) {
						goto F20;
					}
					b1 = b2;
					if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
						return;
					}
					b2 = b2 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
					a2 = a2 * b2;
					if (mc_fmaxf(b2, b1) * hundrd < a2 || cnst1 < a2) {
						goto F20;
					}
				}
F20:
				a2 = a2 * cnst3;
				if (a2 < cnst1) {
					s = gam * (one - mc_sqrtf(a2)) / (a2 + one);
				}
			}
		} else if (dmin == dn2) {
			*ttype = -5;
			 s     = dmin * qurtr;
			 np    = nn - (2 * pp);
			 b1    = mc_blas_vector_at(z, np - 2);
			 b2    = mc_blas_vector_at(z, np - 6);
			 gam   = dn2;

			if (mc_blas_vector_at(z, np - 8) > b2 || mc_blas_vector_at(z, np - 4) > b1) {
				return;
			}
			a2 = mc_blas_vector_at(z, np - 8) / b2 * (mc_blas_vector_at(z, np - 4) / b1 + one);

			if (n0 - i0 > 2) {
				b2 = mc_blas_vector_at(z, nn - 13) / mc_blas_vector_at(z, nn - 15);
				a2 = a2 + b2;
				for (i4 = (nn - 17); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
					if (b2 == zero) {
						goto F40;
					}
					b1 = b2;
					if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
						return;
					}
					b2 = b2 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
					a2 = a2 + b2;
					if (mc_fmaxf(b2, b1) * hundrd < a2 || cnst1 < a2) {
						goto F40;
					}
				}
F40:
				a2 = a2 * cnst3;
			}
			if (a2 < cnst1) {
				s = gam * (one - mc_sqrtf(a2)) / (a2 + one);
			}
		} else {
			if (*ttype == -6) {
				*g = (*g) + ((one - *g) * third);
			} else if (*ttype == -18) {
				*g = qurtr * third;
			} else {
				*g = qurtr;
			}
			 s = *g * dmin;
			*ttype = -6;
		}
	} else if (n0in == n0 + 1) {
		if (dmin1 == dn1 && dmin2 == dn2) {
			*ttype = -7;
			 s     = dmin1 * third;
			if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
				return;
			}
			b1 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
			b2 = b1;
			if (b2 == zero) {
				goto F60;
			}
			for (i4 = ((4 * n0) - 9 + pp); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
				a2 = b1;
				if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
					return;
				}
				b1 = b1 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
				b2 = b2 + b1;
				if (mc_fmaxf(b1, a2) * hundrd < b2) {
					goto F60;
				}
			}
F60:
			b2   = mc_sqrtf(b2 * cnst3);
			a2   = dmin1 / (b2 * b2 + one);
			gap2 = dmin2 * half - a2;
			if (gap2 > zero && gap2 > b2 * a2) {
				s = mc_fmaxf(s, a2 * (one - a2 * cnst2 * (b2 / gap2) * b2));
			} else {
				 s     = mc_fmaxf(s, a2 * (one - b2 * cnst2));
				*ttype = -8;
			}
		} else {
			s = dmin1 * qurtr;
			if (dmin1 == dn1) {
				s = dmin1 * half;
			}
			*ttype = -9;
		}
	} else if (n0in == n0 + 2) {
		if (dmin2 == dn2 && mc_blas_vector_at(z, nn - 5) * two < mc_blas_vector_at(z, nn - 7)) {
			*ttype = -10;
			 s     = dmin2 * third;
			if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
				return;
			}
			b1 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
			b2 = b1;
			if (b2 == zero) {
				goto F80;
			}
			for (i4 = ((4 * n0) - 9 + pp); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
				if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
					return;
				}
				b1 = b1 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
				b2 = b2 + b1;
				if (b1 * hundrd < b2) {
					goto F80;
				}
			}
F80:
			b2   = mc_sqrtf(b2 * cnst3);
			a2   = dmin2 / (b2 * b2 + one);
			gap2 = mc_blas_vector_at(z, nn - 7) + mc_blas_vector_at(z, nn - 9) - mc_sqrtf(mc_blas_vector_at(z, nn - 11)) * mc_sqrtf(mc_blas_vector_at(z, nn - 9)) - a2;
			if (gap2 > zero && gap2 > b2 * a2) {
				s = mc_fmaxf(s, a2 * (one - a2 * cnst2 * (b2 / gap2) * b2));
			} else {
				s = mc_fmaxf(s, a2 * (one - b2 * cnst2));
			}
		} else {
			 s     = dmin2 * qurtr;
			*ttype = -11;
		}
	} else if (n0in > n0 + 2) {
		 s     = zero;
		*ttype = -12;
	}
	*tau = s;
}

#pragma mark - mc_lapack_dlasq4 -

MC_TARGET_PROC void mc_lapack_dlasq4(const int i0, const int n0, const double * z
	, const int pp
	, const int n0in
	, double    dmin
	, double    dmin1
	, double    dmin2
	, double    dn
	, double    dn1
	, double    dn2
	, double *  tau
	, int *     ttype
	, double *  g
) {
	const double hundrd = 100.0, two = 2.0, one = 1.0, half = 0.5, third = 0.333, qurtr = 0.25, zero = 0.0;
	const double cnst1 = 0.563, cnst2 = 1.01, cnst3 = 1.05;

	int i4, nn, np;
	double a2, b1, b2, gam, gap1, gap2, s;

	if (dmin <= zero) {
		*tau   = -(dmin);
		*ttype = -1;
		return;
	}

	s  = zero;
	nn = (4 * n0) + pp;
	if (n0in == n0) {
		if (dmin == dn || dmin == dn1) {
			b1 = mc_sqrt(mc_blas_vector_at(z, nn - 3)) * mc_sqrt(mc_blas_vector_at(z, nn - 5));
			b2 = mc_sqrt(mc_blas_vector_at(z, nn - 7)) * mc_sqrt(mc_blas_vector_at(z, nn - 9));
			a2 = mc_blas_vector_at(z, nn - 7) + mc_blas_vector_at(z, nn - 5);
			if (dmin == dn && dmin1 == dn1) {
				gap2 = dmin2 - a2 - dmin2 * qurtr;
				if (gap2 > zero && gap2 > b2) {
					gap1 = a2 - dn - b2 / gap2 * b2;
				} else {
					gap1 = a2 - dn - (b1 + b2);
				}
				if (gap1 > zero && gap1 > b1) {
					 s     = mc_fmax(dn - b1 / gap1 * b1, dmin * half);
					*ttype = -2;
				} else {
					s = zero;
					if (dn > b1) {
						s = dn - b1;
					}
					if (a2 > b1 + b2) {
						s = mc_fmin(s, a2 - (b1 + b2));
					}
					 s     = mc_fmax(s, dmin * third);
					*ttype = -3;
				}
			} else {
				*ttype = -4;
				s = dmin * qurtr;
				if (dmin == dn) {
					gam = dn;
					a2  = zero;
					if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
						return;
					}
					b2 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
					np = nn - 9;
				} else {
					np  = nn - (2 * pp);
					gam = dn1;
					if (mc_blas_vector_at(z, np - 4) > mc_blas_vector_at(z, np - 2)) {
						return;
					}
					a2 = mc_blas_vector_at(z, np - 4) / mc_blas_vector_at(z, np - 2);
					if (mc_blas_vector_at(z, nn - 9) > mc_blas_vector_at(z, nn - 11)) {
						return;
					}
					b2 = mc_blas_vector_at(z, nn - 9) / mc_blas_vector_at(z, nn - 11);
					np = nn - 13;
				}
				a2 = a2 + b2;
				for (i4 = np; i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
					if (b2 == zero) {
						goto F20;
					}
					b1 = b2;
					if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
						return;
					}
					b2 = b2 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
					a2 = a2 * b2;
					if (mc_fmax(b2, b1) * hundrd < a2 || cnst1 < a2) {
						goto F20;
					}
				}
F20:
				a2 = a2 * cnst3;
				if (a2 < cnst1) {
					s = gam * (one - mc_sqrt(a2)) / (a2 + one);
				}
			}
		} else if (dmin == dn2) {
			*ttype = -5;
			 s     = dmin * qurtr;
			 np    = nn - (2 * pp);
			 b1    = mc_blas_vector_at(z, np - 2);
			 b2    = mc_blas_vector_at(z, np - 6);
			 gam   = dn2;

			if (mc_blas_vector_at(z, np - 8) > b2 || mc_blas_vector_at(z, np - 4) > b1) {
				return;
			}
			a2 = mc_blas_vector_at(z, np - 8) / b2 * (mc_blas_vector_at(z, np - 4) / b1 + one);

			if (n0 - i0 > 2) {
				b2 = mc_blas_vector_at(z, nn - 13) / mc_blas_vector_at(z, nn - 15);
				a2 = a2 + b2;
				for (i4 = (nn - 17); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
					if (b2 == zero) {
						goto F40;
					}
					b1 = b2;
					if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
						return;
					}
					b2 = b2 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
					a2 = a2 + b2;
					if (mc_fmax(b2, b1) * hundrd < a2 || cnst1 < a2) {
						goto F40;
					}
				}
F40:
				a2 = a2 * cnst3;
			}
			if (a2 < cnst1) {
				s = gam * (one - mc_sqrt(a2)) / (a2 + one);
			}
		} else {
			if (*ttype == -6) {
				*g = (*g) + ((one - *g) * third);
			} else if (*ttype == -18) {
				*g = qurtr * third;
			} else {
				*g = qurtr;
			}
			 s = *g * dmin;
			*ttype = -6;
		}
	} else if (n0in == n0 + 1) {
		if (dmin1 == dn1 && dmin2 == dn2) {
			*ttype = -7;
			 s     = dmin1 * third;
			if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
				return;
			}
			b1 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
			b2 = b1;
			if (b2 == zero) {
				goto F60;
			}
			for (i4 = ((4 * n0) - 9 + pp); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
				a2 = b1;
				if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
					return;
				}
				b1 = b1 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
				b2 = b2 + b1;
				if (mc_fmax(b1, a2) * hundrd < b2) {
					goto F60;
				}
			}
F60:
			b2   = mc_sqrt(b2 * cnst3);
			a2   = dmin1 / (b2 * b2 + one);
			gap2 = dmin2 * half - a2;
			if (gap2 > zero && gap2 > b2 * a2) {
				s = mc_fmax(s, a2 * (one - a2 * cnst2 * (b2 / gap2) * b2));
			} else {
				 s     = mc_fmax(s, a2 * (one - b2 * cnst2));
				*ttype = -8;
			}
		} else {
			s = dmin1 * qurtr;
			if (dmin1 == dn1) {
				s = dmin1 * half;
			}
			*ttype = -9;
		}
	} else if (n0in == n0 + 2) {
		if (dmin2 == dn2 && mc_blas_vector_at(z, nn - 5) * two < mc_blas_vector_at(z, nn - 7)) {
			*ttype = -10;
			 s     = dmin2 * third;
			if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
				return;
			}
			b1 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
			b2 = b1;
			if (b2 == zero) {
				goto F80;
			}
			for (i4 = ((4 * n0) - 9 + pp); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
				if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
					return;
				}
				b1 = b1 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
				b2 = b2 + b1;
				if (b1 * hundrd < b2) {
					goto F80;
				}
			}
F80:
			b2   = mc_sqrt(b2 * cnst3);
			a2   = dmin2 / (b2 * b2 + one);
			gap2 = mc_blas_vector_at(z, nn - 7) + mc_blas_vector_at(z, nn - 9) - mc_sqrt(mc_blas_vector_at(z, nn - 11)) * mc_sqrt(mc_blas_vector_at(z, nn - 9)) - a2;
			if (gap2 > zero && gap2 > b2 * a2) {
				s = mc_fmax(s, a2 * (one - a2 * cnst2 * (b2 / gap2) * b2));
			} else {
				s = mc_fmax(s, a2 * (one - b2 * cnst2));
			}
		} else {
			 s     = dmin2 * qurtr;
			*ttype = -11;
		}
	} else if (n0in > n0 + 2) {
		 s     = zero;
		*ttype = -12;
	}
	*tau = s;
}

#pragma mark - mc_lapack_llasq4 -

MC_TARGET_PROC void mc_lapack_llasq4(const int i0, const int n0, const long double * z
	, const int     pp
	, const int     n0in
	, long double   dmin
	, long double   dmin1
	, long double   dmin2
	, long double   dn
	, long double   dn1
	, long double   dn2
	, long double * tau
	, int *         ttype
	, long double * g
) {
	const long double hundrd = 100.0L, two = 2.0L, one = 1.0L, half = 0.5L, third = 0.333L, qurtr = 0.25L, zero = 0.0L;
	const long double cnst1 = 0.563L, cnst2 = 1.01L, cnst3 = 1.05L;

	int i4, nn, np;
	long double a2, b1, b2, gam, gap1, gap2, s;

	if (dmin <= zero) {
		*tau   = -(dmin);
		*ttype = -1;
		return;
	}

	s  = zero;
	nn = (4 * n0) + pp;
	if (n0in == n0) {
		if (dmin == dn || dmin == dn1) {
			b1 = mc_sqrtl(mc_blas_vector_at(z, nn - 3)) * mc_sqrtl(mc_blas_vector_at(z, nn - 5));
			b2 = mc_sqrtl(mc_blas_vector_at(z, nn - 7)) * mc_sqrtl(mc_blas_vector_at(z, nn - 9));
			a2 = mc_blas_vector_at(z, nn - 7) + mc_blas_vector_at(z, nn - 5);
			if (dmin == dn && dmin1 == dn1) {
				gap2 = dmin2 - a2 - dmin2 * qurtr;
				if (gap2 > zero && gap2 > b2) {
					gap1 = a2 - dn - b2 / gap2 * b2;
				} else {
					gap1 = a2 - dn - (b1 + b2);
				}
				if (gap1 > zero && gap1 > b1) {
					 s     = mc_fmaxl(dn - b1 / gap1 * b1, dmin * half);
					*ttype = -2;
				} else {
					s = zero;
					if (dn > b1) {
						s = dn - b1;
					}
					if (a2 > b1 + b2) {
						s = mc_fminl(s, a2 - (b1 + b2));
					}
					 s     = mc_fmaxl(s, dmin * third);
					*ttype = -3;
				}
			} else {
				*ttype = -4;
				s = dmin * qurtr;
				if (dmin == dn) {
					gam = dn;
					a2  = zero;
					if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
						return;
					}
					b2 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
					np = nn - 9;
				} else {
					np  = nn - (2 * pp);
					gam = dn1;
					if (mc_blas_vector_at(z, np - 4) > mc_blas_vector_at(z, np - 2)) {
						return;
					}
					a2 = mc_blas_vector_at(z, np - 4) / mc_blas_vector_at(z, np - 2);
					if (mc_blas_vector_at(z, nn - 9) > mc_blas_vector_at(z, nn - 11)) {
						return;
					}
					b2 = mc_blas_vector_at(z, nn - 9) / mc_blas_vector_at(z, nn - 11);
					np = nn - 13;
				}
				a2 = a2 + b2;
				for (i4 = np; i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
					if (b2 == zero) {
						goto F20;
					}
					b1 = b2;
					if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
						return;
					}
					b2 = b2 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
					a2 = a2 * b2;
					if (mc_fmaxl(b2, b1) * hundrd < a2 || cnst1 < a2) {
						goto F20;
					}
				}
F20:
				a2 = a2 * cnst3;
				if (a2 < cnst1) {
					s = gam * (one - mc_sqrtl(a2)) / (a2 + one);
				}
			}
		} else if (dmin == dn2) {
			*ttype = -5;
			 s     = dmin * qurtr;
			 np    = nn - (2 * pp);
			 b1    = mc_blas_vector_at(z, np - 2);
			 b2    = mc_blas_vector_at(z, np - 6);
			 gam   = dn2;

			if (mc_blas_vector_at(z, np - 8) > b2 || mc_blas_vector_at(z, np - 4) > b1) {
				return;
			}
			a2 = mc_blas_vector_at(z, np - 8) / b2 * (mc_blas_vector_at(z, np - 4) / b1 + one);

			if (n0 - i0 > 2) {
				b2 = mc_blas_vector_at(z, nn - 13) / mc_blas_vector_at(z, nn - 15);
				a2 = a2 + b2;
				for (i4 = (nn - 17); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
					if (b2 == zero) {
						goto F40;
					}
					b1 = b2;
					if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
						return;
					}
					b2 = b2 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
					a2 = a2 + b2;
					if (mc_fmaxl(b2, b1) * hundrd < a2 || cnst1 < a2) {
						goto F40;
					}
				}
F40:
				a2 = a2 * cnst3;
			}
			if (a2 < cnst1) {
				s = gam * (one - mc_sqrtl(a2)) / (a2 + one);
			}
		} else {
			if (*ttype == -6) {
				*g = (*g) + ((one - *g) * third);
			} else if (*ttype == -18) {
				*g = qurtr * third;
			} else {
				*g = qurtr;
			}
			 s = *g * dmin;
			*ttype = -6;
		}
	} else if (n0in == n0 + 1) {
		if (dmin1 == dn1 && dmin2 == dn2) {
			*ttype = -7;
			 s     = dmin1 * third;
			if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
				return;
			}
			b1 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
			b2 = b1;
			if (b2 == zero) {
				goto F60;
			}
			for (i4 = ((4 * n0) - 9 + pp); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
				a2 = b1;
				if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
					return;
				}
				b1 = b1 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
				b2 = b2 + b1;
				if (mc_fmaxl(b1, a2) * hundrd < b2) {
					goto F60;
				}
			}
F60:
			b2   = mc_sqrtl(b2 * cnst3);
			a2   = dmin1 / (b2 * b2 + one);
			gap2 = dmin2 * half - a2;
			if (gap2 > zero && gap2 > b2 * a2) {
				s = mc_fmaxl(s, a2 * (one - a2 * cnst2 * (b2 / gap2) * b2));
			} else {
				 s     = mc_fmaxl(s, a2 * (one - b2 * cnst2));
				*ttype = -8;
			}
		} else {
			s = dmin1 * qurtr;
			if (dmin1 == dn1) {
				s = dmin1 * half;
			}
			*ttype = -9;
		}
	} else if (n0in == n0 + 2) {
		if (dmin2 == dn2 && mc_blas_vector_at(z, nn - 5) * two < mc_blas_vector_at(z, nn - 7)) {
			*ttype = -10;
			 s     = dmin2 * third;
			if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 7)) {
				return;
			}
			b1 = mc_blas_vector_at(z, nn - 5) / mc_blas_vector_at(z, nn - 7);
			b2 = b1;
			if (b2 == zero) {
				goto F80;
			}
			for (i4 = ((4 * n0) - 9 + pp); i4 >= ((4 * i0) - 1 + pp); i4 += -4) {
				if (mc_blas_vector_at(z, i4) > mc_blas_vector_at(z, i4 - 2)) {
					return;
				}
				b1 = b1 * (mc_blas_vector_at(z, i4) / mc_blas_vector_at(z, i4 - 2));
				b2 = b2 + b1;
				if (b1 * hundrd < b2) {
					goto F80;
				}
			}
F80:
			b2   = mc_sqrtl(b2 * cnst3);
			a2   = dmin2 / (b2 * b2 + one);
			gap2 = mc_blas_vector_at(z, nn - 7) + mc_blas_vector_at(z, nn - 9) - mc_sqrtl(mc_blas_vector_at(z, nn - 11)) * mc_sqrtl(mc_blas_vector_at(z, nn - 9)) - a2;
			if (gap2 > zero && gap2 > b2 * a2) {
				s = mc_fmaxl(s, a2 * (one - a2 * cnst2 * (b2 / gap2) * b2));
			} else {
				s = mc_fmaxl(s, a2 * (one - b2 * cnst2));
			}
		} else {
			 s     = dmin2 * qurtr;
			*ttype = -11;
		}
	} else if (n0in > n0 + 2) {
		 s     = zero;
		*ttype = -12;
	}
	*tau = s;
}

#endif /* !MC_LAPACKE_LASQ4_H */

/* EOF */