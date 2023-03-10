//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_minmaxmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MINMAXMX1_H
#define MC_MINMAXMX1_H

#pragma mark - mc_minmaxmx1 -

MC_TARGET_FUNC void mc_minmaxmx1f(const int m, const int n, const int j, const float * a, float * min, float * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p);
	const int wantq   = mc_nonnullptr(q);

	int i = 2;

	if (m > 0) {
		if (m == 1) {
			if (wantmin) {
				*min = a[j];
			}
			if (wantmax) {
				*max = a[j];
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (a[j] > a[n + j]) {
				if (wantmin) {
					*min = a[n + j];
				}
				if (wantmax) {
					*max = a[j];
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = a[j];
				}
				if (wantmax) {
					*max = a[n + j];
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < m; i++) {
				if (a[(n * i) + j] < *min) {
					if (wantmin) {
						*min = a[(n * i) + j];
					}
					if (wantp) {
						*p = i;
					}
				}
				if (a[(n * i) + j] > *max) {
					if (wantmax) {
						*max = a[(n * i) + j];
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_minmaxmx1ff(const int m, const int n, const int j, const float * a, double * min, double * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p);
	const int wantq   = mc_nonnullptr(q);

	int i = 2;

	if (m > 0) {
		if (m == 1) {
			if (wantmin) {
				*min = mc_cast(double, a[j]);
			}
			if (wantmax) {
				*max = mc_cast(double, a[j]);
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (a[j] > a[n + j]) {
				if (wantmin) {
					*min = mc_cast(double, a[n + j]);
				}
				if (wantmax) {
					*max = mc_cast(double, a[j]);
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = mc_cast(double, a[j]);
				}
				if (wantmax) {
					*max = mc_cast(double, a[n + j]);
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < m; i++) {
				if (a[(n * i) + j] < *min) {
					if (wantmin) {
						*min = mc_cast(double, a[(n * i) + j]);
					}
					if (wantp) {
						*p = i;
					}
				}
				if (a[(n * i) + j] > *max) {
					if (wantmax) {
						*max = mc_cast(double, a[(n * i) + j]);
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_minmaxmx1(const int m, const int n, const int j, const double * a, double * min, double * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p);
	const int wantq   = mc_nonnullptr(q);

	int i = 2;

	if (m > 0) {
		if (m == 1) {
			if (wantmin) {
				*min = a[j];
			}
			if (wantmax) {
				*max = a[j];
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (a[j] > a[n + j]) {
				if (wantmin) {
					*min = a[n + j];
				}
				if (wantmax) {
					*max = a[j];
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = a[j];
				}
				if (wantmax) {
					*max = a[n + j];
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < m; i++) {
				if (a[(n * i) + j] < *min) {
					if (wantmin) {
						*min = a[(n * i) + j];
					}
					if (wantp) {
						*p = i;
					}
				}
				if (a[(n * i) + j] > *max) {
					if (wantmax) {
						*max = a[(n * i) + j];
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_minmaxmx1l(const int m, const int n, const int j, const long double * a, long double * min, long double * max, int * p, int * q)
{
	const int wantmin = mc_nonnullptr(min);
	const int wantmax = mc_nonnullptr(max);
	const int wantp   = mc_nonnullptr(p);
	const int wantq   = mc_nonnullptr(q);

	int i = 2;

	if (m > 0) {
		if (m == 1) {
			if (wantmin) {
				*min = a[j];
			}
			if (wantmax) {
				*max = a[j];
			}
			if (wantp) {
				*p = 0;
			}
			if (wantq) {
				*q = 0;
			}
		} else {
			if (a[j] > a[n + j]) {
				if (wantmin) {
					*min = a[n + j];
				}
				if (wantmax) {
					*max = a[j];
				}
				if (wantp) {
					*p = 1;
				}
				if (wantq) {
					*q = 0;
				}
			} else {
				if (wantmin) {
					*min = a[j];
				}
				if (wantmax) {
					*max = a[n + j];
				}
				if (wantp) {
					*p = 0;
				}
				if (wantq) {
					*q = 1;
				}
			}
			for (; i < m; i++) {
				if (a[(n * i) + j] < *min) {
					if (wantmin) {
						*min = a[(n * i) + j];
					}
					if (wantp) {
						*p = i;
					}
				}
				if (a[(n * i) + j] > *max) {
					if (wantmax) {
						*max = a[(n * i) + j];
					}
					if (wantq) {
						*q = i;
					}
				}
			}
		}
	}
}

#endif /* !MC_MINMAXMX1_H */

/* EOF */