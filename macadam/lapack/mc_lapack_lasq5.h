//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasq5.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>
#include <macadam/details/math/mc_fmin.h>

#ifndef MC_LAPACKE_LASQ5_H
#define MC_LAPACKE_LASQ5_H

#pragma mark - mc_lapack_slasq5 -

MC_TARGET_PROC void mc_lapack_slasq5(const int i0, const int n0, float * z, const int pp, float tau, float sigma
	, float *   dmin
	, float *   dmin1
	, float *   dmin2
	, float *   dn
	, float *   dnm1
	, float *   dnm2
	, const int ieee
	, float     eps
) {
	const float half = 0.5f, zero = 0.0f;

	int j4, j4p2;
	float d, emin, temp, dthresh;

	if ((n0 - i0 - 1) <= 0) {
		return;
	}

	dthresh = eps * (sigma + tau);
	if (tau < dthresh * half) {
		tau = zero;
	}
	if (tau != zero) {
		 j4    = (4 * i0) + pp - 3;
		 emin  = mc_blas_vector_at(z, j4 + 4);
		 d     = mc_blas_vector_at(z, j4) - tau;
		*dmin  = d;
		*dmin1 = -mc_blas_vector_at(z, j4);

		if (ieee) {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					 temp                        = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
					 d                           = d * temp - tau;
					*dmin                        = mc_fminf(*dmin, d);
					mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4 - 1) * temp;
					 emin                        = mc_fminf(mc_blas_vector_at(z, j4), emin);
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					 temp                        = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
					 d                           = d * temp - tau;
					*dmin                        = mc_fminf(*dmin, d);
					mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
					 emin                        = mc_fminf(mc_blas_vector_at(z, j4 - 1), emin);
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dnm1                        = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminf(*dmin, *dnm1);

			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dn                          = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminf(*dmin, *dn);
		} else {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
						d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2)) - tau;
					}
					*dmin = mc_fminf(*dmin,d);
					 emin = mc_fminf(emin, mc_blas_vector_at(z, j4));
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
						d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3)) - tau;
					}
					*dmin = mc_fminf(*dmin, d);
					 emin = mc_fminf(emin, mc_blas_vector_at(z, j4 - 1));
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			if (*dnm2 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}

			*dmin                        = mc_fminf(*dmin, *dnm1);
			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			if (*dnm1 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}
			*dmin = mc_fminf(*dmin,*dn);
		}
	} else {
		 j4    = (4 * i0) + pp - 3;
		 emin  = mc_blas_vector_at(z, j4 + 4);
		 d     = mc_blas_vector_at(z, j4) - tau;
		*dmin  = d;
		*dmin1 = -mc_blas_vector_at(z, j4);
		if (ieee) {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					 temp                        = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
					 d                           = d * temp - tau;
					if (d < dthresh) {
						d = zero;
					}
					*dmin                        = mc_fminf(*dmin, d);
					mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4 - 1) * temp;
					 emin                        = mc_fminf(mc_blas_vector_at(z, j4), emin);
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					 temp                        = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
					 d = d * temp - tau;
					if (d < dthresh) {
						d = zero;
					}
					*dmin                        = mc_fminf(*dmin, d);
					mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
					 emin                        = mc_fminf(mc_blas_vector_at(z, j4 - 1), emin);
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dnm1                        = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminf(*dmin, *dnm1);

			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dn                          = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminf(*dmin, *dn);
		} else {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
						d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2)) - tau;
					}
					if (d < dthresh) {
						d = zero;
					}
					*dmin = mc_fminf(*dmin, d);
					 emin = mc_fminf(emin, mc_blas_vector_at(z, j4));
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
						d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3)) - tau;
					}
					if (d < dthresh) {
						d = zero;
					}
					*dmin = mc_fminf(*dmin,d);
					 emin = mc_fminf(emin, mc_blas_vector_at(z, j4 - 1));
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			if (*dnm2 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}

			*dmin                        = mc_fminf(*dmin,*dnm1);
			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			if (*dnm1 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}
			*dmin = mc_fminf(*dmin,*dn);
		}
	}
	mc_blas_vector_at(z, j4 + 2)        = *dn;
	mc_blas_vector_at(z, (4 * n0) - pp) =  emin;
}

#pragma mark - mc_lapack_dlasq5 -

MC_TARGET_PROC void mc_lapack_dlasq5(const int i0, const int n0, double * z, const int pp, double tau, double sigma
	, double *  dmin
	, double *  dmin1
	, double *  dmin2
	, double *  dn
	, double *  dnm1
	, double *  dnm2
	, const int ieee
	, double    eps
) {
	const double half = 0.5, zero = 0.0;

	int j4, j4p2;
	double d, emin, temp, dthresh;

	if ((n0 - i0 - 1) <= 0) {
		return;
	}

	dthresh = eps * (sigma + tau);
	if (tau < dthresh * half) {
		tau = zero;
	}
	if (tau != zero) {
		 j4    = (4 * i0) + pp - 3;
		 emin  = mc_blas_vector_at(z, j4 + 4);
		 d     = mc_blas_vector_at(z, j4) - tau;
		*dmin  = d;
		*dmin1 = -mc_blas_vector_at(z, j4);

		if (ieee) {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					 temp                        = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
					 d                           = d * temp - tau;
					*dmin                        = mc_fmin(*dmin, d);
					mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4 - 1) * temp;
					 emin                        = mc_fmin(mc_blas_vector_at(z, j4), emin);
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					 temp                        = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
					 d                           = d * temp - tau;
					*dmin                        = mc_fmin(*dmin, d);
					mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
					 emin                        = mc_fmin(mc_blas_vector_at(z, j4 - 1), emin);
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dnm1                        = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fmin(*dmin, *dnm1);

			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dn                          = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fmin(*dmin, *dn);
		} else {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
						d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2)) - tau;
					}
					*dmin = mc_fmin(*dmin,d);
					 emin = mc_fmin(emin, mc_blas_vector_at(z, j4));
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
						d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3)) - tau;
					}
					*dmin = mc_fmin(*dmin, d);
					 emin = mc_fmin(emin, mc_blas_vector_at(z, j4 - 1));
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			if (*dnm2 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}

			*dmin                        = mc_fmin(*dmin, *dnm1);
			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			if (*dnm1 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}
			*dmin = mc_fmin(*dmin,*dn);
		}
	} else {
		 j4    = (4 * i0) + pp - 3;
		 emin  = mc_blas_vector_at(z, j4 + 4);
		 d     = mc_blas_vector_at(z, j4) - tau;
		*dmin  = d;
		*dmin1 = -mc_blas_vector_at(z, j4);
		if (ieee) {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					 temp                        = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
					 d                           = d * temp - tau;
					if (d < dthresh) {
						d = zero;
					}
					*dmin                        = mc_fmin(*dmin, d);
					mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4 - 1) * temp;
					 emin                        = mc_fmin(mc_blas_vector_at(z, j4), emin);
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					 temp                        = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
					 d = d * temp - tau;
					if (d < dthresh) {
						d = zero;
					}
					*dmin                        = mc_fmin(*dmin, d);
					mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
					 emin                        = mc_fmin(mc_blas_vector_at(z, j4 - 1), emin);
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dnm1                        = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fmin(*dmin, *dnm1);

			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dn                          = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fmin(*dmin, *dn);
		} else {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
						d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2)) - tau;
					}
					if (d < dthresh) {
						d = zero;
					}
					*dmin = mc_fmin(*dmin, d);
					 emin = mc_fmin(emin, mc_blas_vector_at(z, j4));
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
						d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3)) - tau;
					}
					if (d < dthresh) {
						d = zero;
					}
					*dmin = mc_fmin(*dmin,d);
					 emin = mc_fmin(emin, mc_blas_vector_at(z, j4 - 1));
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			if (*dnm2 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}

			*dmin                        = mc_fmin(*dmin,*dnm1);
			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			if (*dnm1 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}
			*dmin = mc_fmin(*dmin,*dn);
		}
	}
	mc_blas_vector_at(z, j4 + 2)        = *dn;
	mc_blas_vector_at(z, (4 * n0) - pp) =  emin;
}

#pragma mark - mc_lapack_llasq5 -

MC_TARGET_PROC void mc_lapack_llasq5(const int i0, const int n0, long double * z, const int pp, long double tau, long double sigma
	, long double * dmin
	, long double * dmin1
	, long double * dmin2
	, long double * dn
	, long double * dnm1
	, long double * dnm2
	, const int     ieee
	, long double   eps
) {
	const long double half = 0.5L, zero = 0.0L;

	int j4, j4p2;
	long double d, emin, temp, dthresh;

	if ((n0 - i0 - 1) <= 0) {
		return;
	}

	dthresh = eps * (sigma + tau);
	if (tau < dthresh * half) {
		tau = zero;
	}
	if (tau != zero) {
		 j4    = (4 * i0) + pp - 3;
		 emin  = mc_blas_vector_at(z, j4 + 4);
		 d     = mc_blas_vector_at(z, j4) - tau;
		*dmin  = d;
		*dmin1 = -mc_blas_vector_at(z, j4);

		if (ieee) {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					 temp                        = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
					 d                           = d * temp - tau;
					*dmin                        = mc_fminl(*dmin, d);
					mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4 - 1) * temp;
					 emin                        = mc_fminl(mc_blas_vector_at(z, j4), emin);
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					 temp                        = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
					 d                           = d * temp - tau;
					*dmin                        = mc_fminl(*dmin, d);
					mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
					 emin                        = mc_fminl(mc_blas_vector_at(z, j4 - 1), emin);
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dnm1                        = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminl(*dmin, *dnm1);

			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dn                          = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminl(*dmin, *dn);
		} else {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
						d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2)) - tau;
					}
					*dmin = mc_fminl(*dmin,d);
					 emin = mc_fminl(emin, mc_blas_vector_at(z, j4));
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
						d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3)) - tau;
					}
					*dmin = mc_fminl(*dmin, d);
					 emin = mc_fminl(emin, mc_blas_vector_at(z, j4 - 1));
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			if (*dnm2 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}

			*dmin                        = mc_fminl(*dmin, *dnm1);
			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			if (*dnm1 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}
			*dmin = mc_fminl(*dmin,*dn);
		}
	} else {
		 j4    = (4 * i0) + pp - 3;
		 emin  = mc_blas_vector_at(z, j4 + 4);
		 d     = mc_blas_vector_at(z, j4) - tau;
		*dmin  = d;
		*dmin1 = -mc_blas_vector_at(z, j4);
		if (ieee) {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					 temp                        = mc_blas_vector_at(z, j4 + 1) / mc_blas_vector_at(z, j4 - 2);
					 d                           = d * temp - tau;
					if (d < dthresh) {
						d = zero;
					}
					*dmin                        = mc_fminl(*dmin, d);
					mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4 - 1) * temp;
					 emin                        = mc_fminl(mc_blas_vector_at(z, j4), emin);
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					 temp                        = mc_blas_vector_at(z, j4 + 2) / mc_blas_vector_at(z, j4 - 3);
					 d = d * temp - tau;
					if (d < dthresh) {
						d = zero;
					}
					*dmin                        = mc_fminl(*dmin, d);
					mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4) * temp;
					 emin                        = mc_fminl(mc_blas_vector_at(z, j4 - 1), emin);
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dnm1                        = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminl(*dmin, *dnm1);

			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			mc_blas_vector_at(z, j4)     = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
			*dn                          = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			*dmin                        = mc_fminl(*dmin, *dn);
		} else {
			if (pp == 0) {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 2) = d + mc_blas_vector_at(z, j4 - 1);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4 + 1) * (mc_blas_vector_at(z, j4 - 1) / mc_blas_vector_at(z, j4 - 2));
						d                        = mc_blas_vector_at(z, j4 + 1) * (d / mc_blas_vector_at(z, j4 - 2)) - tau;
					}
					if (d < dthresh) {
						d = zero;
					}
					*dmin = mc_fminl(*dmin, d);
					 emin = mc_fminl(emin, mc_blas_vector_at(z, j4));
				}
			} else {
				for (j4 = (4 * i0); j4 <= (4 * (n0 - 3)); j4 += 4) {
					mc_blas_vector_at(z, j4 - 3) = d + mc_blas_vector_at(z, j4);
					if (d < zero) {
						return;
					} else {
						mc_blas_vector_at(z, j4 - 1) = mc_blas_vector_at(z, j4 + 2) * (mc_blas_vector_at(z, j4) / mc_blas_vector_at(z, j4 - 3));
						d                            = mc_blas_vector_at(z, j4 + 2) * (d / mc_blas_vector_at(z, j4 - 3)) - tau;
					}
					if (d < dthresh) {
						d = zero;
					}
					*dmin = mc_fminl(*dmin,d);
					 emin = mc_fminl(emin, mc_blas_vector_at(z, j4 - 1));
				}
			}

			*dnm2                        =  d;
			*dmin2                       = *dmin;
			 j4                          =  (4 * (n0 - 2)) - pp;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm2 + mc_blas_vector_at(z, j4p2);
			if (*dnm2 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dnm1                    = mc_blas_vector_at(z, j4p2 + 2) * (*dnm2 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}

			*dmin                        = mc_fminl(*dmin,*dnm1);
			*dmin1                       = *dmin;
			 j4                          =  j4 + 4;
			 j4p2                        =  j4 + (2 * pp) - 1;
			mc_blas_vector_at(z, j4 - 2) = *dnm1 + mc_blas_vector_at(z, j4p2);
			if (*dnm1 < zero) {
				return;
			} else {
				mc_blas_vector_at(z, j4) = mc_blas_vector_at(z, j4p2 + 2) * (mc_blas_vector_at(z, j4p2) / mc_blas_vector_at(z, j4 - 2));
				*dn                      = mc_blas_vector_at(z, j4p2 + 2) * (*dnm1 / mc_blas_vector_at(z, j4 - 2)) - tau;
			}
			*dmin = mc_fminl(*dmin,*dn);
		}
	}
	mc_blas_vector_at(z, j4 + 2)        = *dn;
	mc_blas_vector_at(z, (4 * n0) - pp) =  emin;
}

#endif /* !MC_LAPACKE_LASQ5_H */

/* EOF */