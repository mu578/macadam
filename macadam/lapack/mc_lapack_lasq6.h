//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasq6.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/details/math/mc_fmin.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_LAPACKE_LASQ6_H
#define MC_LAPACKE_LASQ6_H

#pragma mark - mc_lapack_slasq6 -

MC_TARGET_PROC void mc_lapack_slasq6(const int i0, const int n0, float * z, const int pp
	, float * dmin
	, float * dmin1
	, float * dmin2
	, float * dn
	, float * dnm1
	, float * dnm2
) {
	const float zero = 0.0f;

	int j4, j4p2;
	float d, emin, safmin, temp;

	if ((n0 - i0 - 1) <= 0) {
		return;
	}

	 safmin = mc_lapack_slamch('S');
	 j4     = (4 * i0) + pp - 3;
	 emin   = mc_blas_vector_at(z, j4 + 4);
	 d      = mc_blas_vector_at(z, j4);
	*dmin   = d;

	if (pp == 0) {
		for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
			mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
			if (mc_blas_vector_at(z, j4 - 2) == zero) {
				mc_blas_vector_at(z, j4) = zero;
				 d                       = mc_blas_vector_at(z, j4 + 1);
				*dmin                    = d;
				 emin                    = zero;
			} else if (
				   safmin * mc_blas_vector_at(z, j4 + 1) < mc_blas_vector_at(z, j4 - 2)
				&& safmin * mc_blas_vector_at(z, j4 - 2) < mc_blas_vector_at(z, j4 + 1)
			) {
				temp                     = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 - 1) * temp;
				d                        = mc_raise2f(temp);
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
				d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2));
			}
			*dmin = mc_fminf(*dmin, d);
			 emin = mc_fminf(emin, mc_blas_vector_at(z, j4));
		}
	} else {
		for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
			mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
			if (mc_blas_vector_at(z, j4 - 3) == zero) {
				mc_blas_vector_at(z, j4 - 1) = zero;
				 d                           = mc_blas_vector_at(z, j4 + 2);
				*dmin                        = d;
				 emin                        = zero;
			} else if (
				   safmin * mc_blas_vector_at(z, j4 + 2) < mc_blas_vector_at(z, j4 - 3)
				&& safmin * mc_blas_vector_at(z, j4 - 3) < mc_blas_vector_at(z, j4 + 2)
			) {
				temp = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
				mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
				d                            = d * temp;
			} else {
				mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
				d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3));
			}
			*dmin = mc_fminf(*dmin, d);
			 emin = mc_fminf(emin, mc_blas_vector_at(z, j4 - 1));
		}
	}

	*dnm2                        = d;
	*dmin2                       = *dmin;
	 j4                          = 4 * (n0 - 2) - pp;
	 j4p2                        = j4 + (2 * pp) - 1;
	mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
	if (mc_blas_vector_at(z, j4 - 2) == zero) {
		mc_blas_vector_at(z, j4) = zero;
		*dnm1                    = mc_blas_vector_at(z, j4p2 + 2);
		*dmin                    = *dnm1;
		 emin                    = zero;
	} else if (
		   safmin * mc_blas_vector_at(z, j4p2 + 2) < mc_blas_vector_at(z, j4 - 2)
		&& safmin * mc_blas_vector_at(z, j4   - 2) < mc_blas_vector_at(z, j4p2 + 2)
	) {
		 temp                    = mc_blas_vector_at(z, j4p2 + 2) / mc_blas_vector_at(z, j4 - 2);
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2) * temp;
		*dnm1                    = *dnm2 * temp;
	} else {
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
		*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2));
	}
	*dmin                        = mc_fminf(*dmin, *dnm1);
	*dmin1                       = *dmin;
	 j4                          = j4 + 4;
	 j4p2                        = j4 + (2 * pp) - 1;
	mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
	if (mc_blas_vector_at(z, j4 - 2) == zero) {
		mc_blas_vector_at(z, j4) = zero;
		*dn                      = mc_blas_vector_at(z, j4p2 + 2);
		*dmin                    = *dn;
		 emin                    = zero;
	} else if (
		   safmin * mc_blas_vector_at(z, j4p2 + 2) < mc_blas_vector_at(z, j4 - 2)
		&& safmin * mc_blas_vector_at(z, j4   - 2) < mc_blas_vector_at(z, j4p2 + 2)
	) {
		 temp                    = mc_blas_vector_at(z, j4p2 + 2) / mc_blas_vector_at(z, j4 - 2);
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2) * temp;
		*dn = *dnm1 * temp;
	} else {
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
		*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2));
	}
	*dmin                                = mc_fminf(*dmin, *dn);
	mc_blas_vector_at(z, j4 + 2)         = *dn;
	mc_blas_vector_at(z, (4 * n0 ) - pp) = emin;
}

#pragma mark - mc_lapack_dlasq6 -

MC_TARGET_PROC void mc_lapack_dlasq6(const int i0, const int n0, double * z, const int pp
	, double * dmin
	, double * dmin1
	, double * dmin2
	, double * dn
	, double * dnm1
	, double * dnm2
) {
	const double zero = 0.0;

	int j4, j4p2;
	double d, emin, safmin, temp;

	if ((n0 - i0 - 1) <= 0) {
		return;
	}

	 safmin = mc_lapack_dlamch('S');
	 j4     = (4 * i0) + pp - 3;
	 emin   = mc_blas_vector_at(z, j4 + 4);
	 d      = mc_blas_vector_at(z, j4);
	*dmin   = d;

	if (pp == 0) {
		for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
			mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
			if (mc_blas_vector_at(z, j4 - 2) == zero) {
				mc_blas_vector_at(z, j4) = zero;
				 d                       = mc_blas_vector_at(z, j4 + 1);
				*dmin                    = d;
				 emin                    = zero;
			} else if (
				   safmin * mc_blas_vector_at(z, j4 + 1) < mc_blas_vector_at(z, j4 - 2)
				&& safmin * mc_blas_vector_at(z, j4 - 2) < mc_blas_vector_at(z, j4 + 1)
			) {
				temp                     = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 - 1) * temp;
				d                        = mc_raise2(temp);
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
				d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2));
			}
			*dmin = mc_fmin(*dmin, d);
			 emin = mc_fmin(emin, mc_blas_vector_at(z, j4));
		}
	} else {
		for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
			mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
			if (mc_blas_vector_at(z, j4 - 3) == zero) {
				mc_blas_vector_at(z, j4 - 1) = zero;
				 d                           = mc_blas_vector_at(z, j4 + 2);
				*dmin                        = d;
				 emin                        = zero;
			} else if (
				   safmin * mc_blas_vector_at(z, j4 + 2) < mc_blas_vector_at(z, j4 - 3)
				&& safmin * mc_blas_vector_at(z, j4 - 3) < mc_blas_vector_at(z, j4 + 2)
			) {
				temp = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
				mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
				d                            = d * temp;
			} else {
				mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
				d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3));
			}
			*dmin = mc_fmin(*dmin, d);
			 emin = mc_fmin(emin, mc_blas_vector_at(z, j4 - 1));
		}
	}

	*dnm2                        = d;
	*dmin2                       = *dmin;
	 j4                          = 4 * (n0 - 2) - pp;
	 j4p2                        = j4 + (2 * pp) - 1;
	mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
	if (mc_blas_vector_at(z, j4 - 2) == zero) {
		mc_blas_vector_at(z, j4) = zero;
		*dnm1                    = mc_blas_vector_at(z, j4p2 + 2);
		*dmin                    = *dnm1;
		 emin                    = zero;
	} else if (
		   safmin * mc_blas_vector_at(z, j4p2 + 2) < mc_blas_vector_at(z, j4 - 2)
		&& safmin * mc_blas_vector_at(z, j4   - 2) < mc_blas_vector_at(z, j4p2 + 2)
	) {
		 temp                    = mc_blas_vector_at(z, j4p2 + 2) / mc_blas_vector_at(z, j4 - 2);
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2) * temp;
		*dnm1                    = *dnm2 * temp;
	} else {
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
		*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2));
	}
	*dmin                        = mc_fmin(*dmin, *dnm1);
	*dmin1                       = *dmin;
	 j4                          = j4 + 4;
	 j4p2                        = j4 + (2 * pp) - 1;
	mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
	if (mc_blas_vector_at(z, j4 - 2) == zero) {
		mc_blas_vector_at(z, j4) = zero;
		*dn                      = mc_blas_vector_at(z, j4p2 + 2);
		*dmin                    = *dn;
		 emin                    = zero;
	} else if (
		   safmin * mc_blas_vector_at(z, j4p2 + 2) < mc_blas_vector_at(z, j4 - 2)
		&& safmin * mc_blas_vector_at(z, j4   - 2) < mc_blas_vector_at(z, j4p2 + 2)
	) {
		 temp                    = mc_blas_vector_at(z, j4p2 + 2) / mc_blas_vector_at(z, j4 - 2);
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2) * temp;
		*dn = *dnm1 * temp;
	} else {
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
		*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2));
	}
	*dmin                                = mc_fmin(*dmin, *dn);
	mc_blas_vector_at(z, j4 + 2)         = *dn;
	mc_blas_vector_at(z, (4 * n0 ) - pp) = emin;
}

#pragma mark - mc_lapack_llasq6 -

MC_TARGET_PROC void mc_lapack_llasq6(const int i0, const int n0, long double * z, const int pp
	, long double * dmin
	, long double * dmin1
	, long double * dmin2
	, long double * dn
	, long double * dnm1
	, long double * dnm2
) {
	const long double zero = 0.0L;

	int j4, j4p2;
	long double d, emin, safmin, temp;

	if ((n0 - i0 - 1) <= 0) {
		return;
	}

	 safmin = mc_lapack_llamch('S');
	 j4     = (4 * i0) + pp - 3;
	 emin   = mc_blas_vector_at(z, j4 + 4);
	 d      = mc_blas_vector_at(z, j4);
	*dmin   = d;

	if (pp == 0) {
		for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
			mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
			if (mc_blas_vector_at(z, j4 - 2) == zero) {
				mc_blas_vector_at(z, j4) = zero;
				 d                       = mc_blas_vector_at(z, j4 + 1);
				*dmin                    = d;
				 emin                    = zero;
			} else if (
				   safmin * mc_blas_vector_at(z, j4 + 1) < mc_blas_vector_at(z, j4 - 2)
				&& safmin * mc_blas_vector_at(z, j4 - 2) < mc_blas_vector_at(z, j4 + 1)
			) {
				temp                     = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 - 1) * temp;
				d                        = mc_raise2l(temp);
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
				d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2));
			}
			*dmin = mc_fminl(*dmin, d);
			 emin = mc_fminl(emin, mc_blas_vector_at(z, j4));
		}
	} else {
		for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
			mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
			if (mc_blas_vector_at(z, j4 - 3) == zero) {
				mc_blas_vector_at(z, j4 - 1) = zero;
				 d                           = mc_blas_vector_at(z, j4 + 2);
				*dmin                        = d;
				 emin                        = zero;
			} else if (
				   safmin * mc_blas_vector_at(z, j4 + 2) < mc_blas_vector_at(z, j4 - 3)
				&& safmin * mc_blas_vector_at(z, j4 - 3) < mc_blas_vector_at(z, j4 + 2)
			) {
				temp = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
				mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
				d                            = d * temp;
			} else {
				mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
				d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3));
			}
			*dmin = mc_fminl(*dmin, d);
			 emin = mc_fminl(emin, mc_blas_vector_at(z, j4 - 1));
		}
	}

	*dnm2                        = d;
	*dmin2                       = *dmin;
	 j4                          = 4 * (n0 - 2) - pp;
	 j4p2                        = j4 + (2 * pp) - 1;
	mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
	if (mc_blas_vector_at(z, j4 - 2) == zero) {
		mc_blas_vector_at(z, j4) = zero;
		*dnm1                    = mc_blas_vector_at(z, j4p2 + 2);
		*dmin                    = *dnm1;
		 emin                    = zero;
	} else if (
		   safmin * mc_blas_vector_at(z, j4p2 + 2) < mc_blas_vector_at(z, j4 - 2)
		&& safmin * mc_blas_vector_at(z, j4   - 2) < mc_blas_vector_at(z, j4p2 + 2)
	) {
		 temp                    = mc_blas_vector_at(z, j4p2 + 2) / mc_blas_vector_at(z, j4 - 2);
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2) * temp;
		*dnm1                    = *dnm2 * temp;
	} else {
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
		*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2));
	}
	*dmin                        = mc_fminl(*dmin, *dnm1);
	*dmin1                       = *dmin;
	 j4                          = j4 + 4;
	 j4p2                        = j4 + (2 * pp) - 1;
	mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
	if (mc_blas_vector_at(z, j4 - 2) == zero) {
		mc_blas_vector_at(z, j4) = zero;
		*dn                      = mc_blas_vector_at(z, j4p2 + 2);
		*dmin                    = *dn;
		 emin                    = zero;
	} else if (
		   safmin * mc_blas_vector_at(z, j4p2 + 2) < mc_blas_vector_at(z, j4 - 2)
		&& safmin * mc_blas_vector_at(z, j4   - 2) < mc_blas_vector_at(z, j4p2 + 2)
	) {
		 temp                    = mc_blas_vector_at(z, j4p2 + 2) / mc_blas_vector_at(z, j4 - 2);
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2) * temp;
		*dn = *dnm1 * temp;
	} else {
		mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
		*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2));
	}
	*dmin                                = mc_fminl(*dmin, *dn);
	mc_blas_vector_at(z, j4 + 2)         = *dn;
	mc_blas_vector_at(z, (4 * n0 ) - pp) = emin;
}

#endif /* !MC_LAPACKE_LASQ6_H */

/* EOF */