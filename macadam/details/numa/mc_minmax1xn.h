//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_minmax1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MINMAX1XN_H
#define MC_MINMAX1XN_H

#pragma mark - mc_minmax1xn -

MC_TARGET_FUNC void mc_minmax1xnf(const int n, const float * x, float * min, float * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p) && wantmin;
	const int wantq   = mc_nonnullptr(q) && wantmax;

	int i = 2;

	if (n > 0) {
		if (n == 1) {
			if (wantmin) {
				*min = x[0];
			}
			if (wantmax) {
				*max = x[0];
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (x[0] > x[1]) {
				if (wantmin) {
					*min = x[1];
				}
				if (wantmax) {
					*max = x[0];
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = x[0];
				}
				if (wantmax) {
					*max = x[1];
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < n; i++) {
				if (x[i] < *min) {
					if (wantmin) {
						*min = x[i];
					}
					if (wantp) {
						*p = i;
					}
				}
				if (x[i] > *max) {
					if (wantmax) {
						*max = x[i];
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_minmax1xnff(const int n, const float * x, double * min, double * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p) && wantmin;
	const int wantq   = mc_nonnullptr(q) && wantmax;

	int i = 2;

	if (n > 0) {
		if (n == 1) {
			if (wantmin) {
				*min = mc_cast(double, x[0]);
			}
			if (wantmax) {
				*max = mc_cast(double, x[0]);
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (x[0] > x[1]) {
				if (wantmin) {
					*min = mc_cast(double, x[1]);
				}
				if (wantmax) {
					*max = mc_cast(double, x[0]);
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = mc_cast(double, x[0]);
				}
				if (wantmax) {
					*max = mc_cast(double, x[1]);
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < n; i++) {
				if (x[i] < *min) {
					if (wantmin) {
						*min = mc_cast(double, x[i]);
					}
					if (wantp) {
						*p = i;
					}
				}
				if (x[i] > *max) {
					if (wantmax) {
						*max = mc_cast(double, x[i]);
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_minmax1xn(const int n, const double * x, double * min, double * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p) && wantmin;
	const int wantq   = mc_nonnullptr(q) && wantmax;

	int i = 2;

	if (n > 0) {
		if (n == 1) {
			if (wantmin) {
				*min = x[0];
			}
			if (wantmax) {
				*max = x[0];
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (x[0] > x[1]) {
				if (wantmin) {
					*min = x[1];
				}
				if (wantmax) {
					*max = x[0];
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = x[0];
				}
				if (wantmax) {
					*max = x[1];
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < n; i++) {
				if (x[i] < *min) {
					if (wantmin) {
						*min = x[i];
					}
					if (wantp) {
						*p = i;
					}
				}
				if (x[i] > *max) {
					if (wantmax) {
						*max = x[i];
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_minmax1xnl(const int n, const long double * x, long double * min, long double * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p) && wantmin;
	const int wantq   = mc_nonnullptr(q) && wantmax;

	int i = 2;

	if (n > 0) {
		if (n == 1) {
			if (wantmin) {
				*min = x[0];
			}
			if (wantmax) {
				*max = x[0];
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (x[0] > x[1]) {
				if (wantmin) {
					*min = x[1];
				}
				if (wantmax) {
					*max = x[0];
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = x[0];
				}
				if (wantmax) {
					*max = x[1];
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < n; i++) {
				if (x[i] < *min) {
					if (wantmin) {
						*min = x[i];
					}
					if (wantp) {
						*p = i;
					}
				}
				if (x[i] > *max) {
					if (wantmax) {
						*max = x[i];
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

#endif /* !MC_MINMAX1XN_H */

/* EOF */