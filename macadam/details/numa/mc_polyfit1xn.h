//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_polyfit1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zeros1xn.h>

#ifndef MC_POLYFIT1XN_H
#define MC_POLYFIT1XN_H

#pragma mark - mc_polyfit1xn -

MC_TARGET_FUNC int mc_polyfit1xnf(const int n, const int d, const float * x, const float * y, float * MC_TARGET_RESTRICT w, float * MC_TARGET_RESTRICT c)
{
//!# Requires x[n], y[n], w[d x 6] and c[d + 1] where 0 < d < n.
//!#     n - Number of query points and fitted values i.e x and y respectively.
//!#     d - Order or degree of polynomial fit specified as a positive integer.
//!#     x - The query points.
//!#     y - Fitted values at query points.
//!#         Specifies the polynomial power of the right-most coefficient in vector c.
//!#     w - Working buffer being of size w[d * 6].
//!#     c - Least-squares fit polynomial coefficients, returned as a vector c[d + 1],
//!#         contains the polynomial coefficients in ascending powers, the highest power
//!#         being the last index. Octave yield results in descending powers order.
//!# \note: based on Nate Domin redux approach polyfit function. Todo: full QR based.

	const int e = d + 1;
	const int f = 6 * d;
	const int g = 2 * e;
	const int h = g + 1;

	float * b;
	float * p;
	float * a;

	float xi, yi, sumpx;

	int i = 0, j, k;

	if (n < 1 || !(n > d)) {
		return -1;
	}
	
	mc_zeros1xnf(f, w);

//!# Assigning partitions.
	b    = w + 0;
	p    = w + e;
	a    = w + h;
	p[0] = n;

//!# Building column-vectors, under/overflow unsafe.
	for (; i < n; i++) {
		xi    = x[i];
		yi    = y[i];
		sumpx = 1.0f;
		for (j = 0; j < e; j++) {
			b[j]  = b[j] + (yi * sumpx);
			sumpx = sumpx * xi;
		}
	}

//!# Sum of power of x, under/overflow unsafe.
	for (i = 0; i < n; i++) {
		xi    = x[i];
		sumpx = xi;
		for (j = 1; j < h; j++) {
			p[j]  = p[j] + sumpx;
			sumpx = sumpx * xi;
		}
	}

//!# Initializing the reduction matrix.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			a[(g * i) + j] = p[i + j];
		}
		a[(g * i) + (i + e)] = 1.0f;
	}

//!# Finding inverse of the left side.
	for (i = 0; i < e; i++) {
		xi = a[(g * i) + i];
		if (xi == 0.0f) {
			return -1;
		}
		xi = 1.0f / xi;
		for (k = 0; k < g; k++) {
			a[(g * i) + k] = a[(g * i) + k] * xi;
		}
		for (j = 0; j < e; j++) {
			if ((j - i) != 0) {
				yi = a[(g * i) + i];
				for (k = 0; k < g; k++) {
					a[(g * i) + k] = a[(g * i) + k] - (yi * a[(g * i) + k]);
				}
			}
		}
	}

//!# Solving, assigning coefficients.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			xi = 0.0f;
			for (k = 0; k < e; k++) {
				xi = xi + (a[(g * i) + (k + e)] * b[k]);
			}
			c[i] = xi;
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_polyfit1xnff(const int n, const int d, const float * x, const float * y, double * MC_TARGET_RESTRICT w, double * MC_TARGET_RESTRICT c)
{
//!# Requires x[n], y[n], w[d x 6] and c[d + 1] where 0 < d < n.
//!#     n - Number of query points and fitted values i.e x and y respectively.
//!#     d - Order or degree of polynomial fit specified as a positive integer.
//!#     x - The query points.
//!#     y - Fitted values at query points.
//!#         Specifies the polynomial power of the right-most coefficient in vector c.
//!#     w - Working buffer being of size w[d * 6].
//!#     c - Least-squares fit polynomial coefficients, returned as a vector c[d + 1],
//!#         contains the polynomial coefficients in ascending powers, the highest power
//!#         being the last index. Octave yield results in descending powers order.
//!# \note: based on Nate Domin redux approach polyfit function. Todo: full QR based.

	const int e = d + 1;
	const int f = 6 * d;
	const int g = 2 * e;
	const int h = g + 1;

	double * b;
	double * p;
	double * a;

	double xi, yi, sumpx;

	int i = 0, j, k;

	if (n < 1 || !(n > d)) {
		return -1;
	}
	
	mc_zeros1xn(f, w);

//!# Assigning partitions.
	b    = w + 0;
	p    = w + e;
	a    = w + h;
	p[0] = n;

//!# Building column-vectors, under/overflow unsafe.
	for (; i < n; i++) {
		xi    = mc_cast(double, x[i]);
		yi    = mc_cast(double, y[i]);
		sumpx = 1.0;
		for (j = 0; j < e; j++) {
			b[j]  = b[j] + (yi * sumpx);
			sumpx = sumpx * xi;
		}
	}

//!# Sum of power of x, under/overflow unsafe.
	for (i = 0; i < n; i++) {
		xi    = mc_cast(double, x[i]);
		sumpx = xi;
		for (j = 1; j < h; j++) {
			p[j]  = p[j] + sumpx;
			sumpx = sumpx * xi;
		}
	}

//!# Initializing the reduction matrix.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			a[(g * i) + j] = p[i + j];
		}
		a[(g * i) + (i + e)] = 1.0;
	}

//!# Finding inverse of the left side.
	for (i = 0; i < e; i++) {
		xi = a[(g * i) + i];
		if (xi == 0.0) {
			return -1;
		}
		xi = 1.0 / xi;
		for (k = 0; k < g; k++) {
			a[(g * i) + k] = a[(g * i) + k] * xi;
		}
		for (j = 0; j < e; j++) {
			if ((j - i) != 0) {
				yi = a[(g * i) + i];
				for (k = 0; k < g; k++) {
					a[(g * i) + k] = a[(g * i) + k] - (yi * a[(g * i) + k]);
				}
			}
		}
	}

//!# Solving, assigning coefficients.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			xi = 0.0;
			for (k = 0; k < e; k++) {
				xi = xi + (a[(g * i) + (k + e)] * b[k]);
			}
			c[i] = xi;
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_polyfit1xn(const int n, const int d, const double * x, const double * y, double * MC_TARGET_RESTRICT w, double * MC_TARGET_RESTRICT c)
{
//!# Requires x[n], y[n], w[d x 6] and c[d + 1] where 0 < d < n.
//!#     n - Number of query points and fitted values i.e x and y respectively.
//!#     d - Order or degree of polynomial fit specified as a positive integer.
//!#     x - The query points.
//!#     y - Fitted values at query points.
//!#         Specifies the polynomial power of the right-most coefficient in vector c.
//!#     w - Working buffer being of size w[d * 6].
//!#     c - Least-squares fit polynomial coefficients, returned as a vector c[d + 1],
//!#         contains the polynomial coefficients in ascending powers, the highest power
//!#         being the last index. Octave yield results in descending powers order.
//!# \note: based on Nate Domin redux approach polyfit function. Todo: full QR based.

	const int e = d + 1;
	const int f = 6 * d;
	const int g = 2 * e;
	const int h = g + 1;

	double * b;
	double * p;
	double * a;

	double xi, yi, sumpx;

	int i = 0, j, k;

	if (n < 1 || !(n > d)) {
		return -1;
	}
	
	mc_zeros1xn(f, w);

//!# Assigning partitions.
	b    = w + 0;
	p    = w + e;
	a    = w + h;
	p[0] = n;

//!# Building column-vectors, under/overflow unsafe.
	for (; i < n; i++) {
		xi    = x[i];
		yi    = y[i];
		sumpx = 1.0;
		for (j = 0; j < e; j++) {
			b[j]  = b[j] + (yi * sumpx);
			sumpx = sumpx * xi;
		}
	}

//!# Sum of power of x, under/overflow unsafe.
	for (i = 0; i < n; i++) {
		xi    = x[i];
		sumpx = xi;
		for (j = 1; j < h; j++) {
			p[j]  = p[j] + sumpx;
			sumpx = sumpx * xi;
		}
	}

//!# Initializing the reduction matrix.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			a[(g * i) + j] = p[i + j];
		}
		a[(g * i) + (i + e)] = 1.0;
	}

//!# Finding inverse of the left side.
	for (i = 0; i < e; i++) {
		xi = a[(g * i) + i];
		if (xi == 0.0) {
			return -1;
		}
		xi = 1.0 / xi;
		for (k = 0; k < g; k++) {
			a[(g * i) + k] = a[(g * i) + k] * xi;
		}
		for (j = 0; j < e; j++) {
			if ((j - i) != 0) {
				yi = a[(g * i) + i];
				for (k = 0; k < g; k++) {
					a[(g * i) + k] = a[(g * i) + k] - (yi * a[(g * i) + k]);
				}
			}
		}
	}

//!# Solving, assigning coefficients.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			xi = 0.0;
			for (k = 0; k < e; k++) {
				xi = xi + (a[(g * i) + (k + e)] * b[k]);
			}
			c[i] = xi;
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_polyfit1xnl(const int n, const int d, const long double * x, const long double * y, long double * MC_TARGET_RESTRICT w, long double * MC_TARGET_RESTRICT c)
{
//!# Requires x[n], y[n], w[d x 6] and c[d + 1] where 0 < d < n.
//!#     n - Number of query points and fitted values i.e x and y respectively.
//!#     d - Order or degree of polynomial fit specified as a positive integer.
//!#     x - The query points.
//!#     y - Fitted values at query points.
//!#         Specifies the polynomial power of the right-most coefficient in vector c.
//!#     w - Working buffer being of size w[d * 6].
//!#     c - Least-squares fit polynomial coefficients, returned as a vector c[d + 1],
//!#         contains the polynomial coefficients in ascending powers, the highest power
//!#         being the last index. Octave yield results in descending powers order.
//!# \note: based on Nate Domin redux approach polyfit function. Todo: full QR based.

	const int e = d + 1;
	const int f = 6 * d;
	const int g = 2 * e;
	const int h = g + 1;

	long double * b;
	long double * p;
	long double * a;

	long double xi, yi, sumpx;

	int i = 0, j, k;

	if (n < 1 || !(n > d)) {
		return -1;
	}
	
	mc_zeros1xnl(f, w);

//!# Assigning partitions.
	b    = w + 0;
	p    = w + e;
	a    = w + h;
	p[0] = n;

//!# Building column-vectors, under/overflow unsafe.
	for (; i < n; i++) {
		xi    = x[i];
		yi    = y[i];
		sumpx = 1.0L;
		for (j = 0; j < e; j++) {
			b[j]  = b[j] + (yi * sumpx);
			sumpx = sumpx * xi;
		}
	}

//!# Sum of power of x, under/overflow unsafe.
	for (i = 0; i < n; i++) {
		xi    = x[i];
		sumpx = xi;
		for (j = 1; j < h; j++) {
			p[j]  = p[j] + sumpx;
			sumpx = sumpx * xi;
		}
	}

//!# Initializing the reduction matrix.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			a[(g * i) + j] = p[i + j];
		}
		a[(g * i) + (i + e)] = 1.0L;
	}

//!# Finding inverse of the left side.
	for (i = 0; i < e; i++) {
		xi = a[(g * i) + i];
		if (xi == 0.0L) {
			return -1;
		}
		xi = 1.0L / xi;
		for (k = 0; k < g; k++) {
			a[(g * i) + k] = a[(g * i) + k] * xi;
		}
		for (j = 0; j < e; j++) {
			if ((j - i) != 0) {
				yi = a[(g * i) + i];
				for (k = 0; k < g; k++) {
					a[(g * i) + k] = a[(g * i) + k] - (yi * a[(g * i) + k]);
				}
			}
		}
	}

//!# Solving, assigning coefficients.
	for (i = 0; i < e; i++) {
		for (j = 0; j < e; j++) {
			xi = 0.0L;
			for (k = 0; k < e; k++) {
				xi = xi + (a[(g * i) + (k + e)] * b[k]);
			}
			c[i] = xi;
		}
	}
	return 0;
}

#endif /* !MC_POLYFIT1XN_H */

/* EOF */