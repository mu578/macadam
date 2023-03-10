//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_magicnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MAGICNXN_H
#define MC_MAGICNXN_H

#pragma mark - mc_magicnxn -

MC_TARGET_FUNC int mc_magicnxnf(const int n, float * a)
{
	int i, j, k, ki, kj;
	int r, l;

	float w;

	if (n > 2 && ((n % 2) == 1)) {
//!# Odd order.
		mc_zerosnxnf(n, a);
		i              = 0;
		j              = n / 2;
		k              = 2;
		ki             = i;
		kj             = j;
		a[(n * i) + j] = 1.0f;
		for (; k <= n * n; k++) {
			i = i - 1 < 0  ? n - 1 : i - 1;
			j = j + 1 == n ? 0     : j + 1;
			if (a[(n * i) + j] != 0.0f) {
				i = ki;
				j = kj;
				i = i + 1 == n ? 0 : i + 1;
			}
			a[(n * i) + j] = mc_cast(float, k);
			ki             = i;
			kj             = j; 
		}
		return 0;
	} else if (n > 2 && ((n % 2 == 0 && (n /2 ) % 2 == 1))) { 
//!# Single-even order.
		r = n / 2;
		i = 0;
		j = r / 2;

//!# Odd order partitions.
		for (k = 1; k <= r * r; k++) {
			a[(n * i) + j]             = mc_cast(float, k            );
			a[(n * (i + r)) + (j + r)] = mc_cast(float, k +     r * r);
			a[(n * i) + (j + r)]       = mc_cast(float, k + 2 * r * r);
			a[(n * (i + r)) + j]       = mc_cast(float, k + 3 * r * r);
			if (!(k % r == 0)) {
				i = (i == 0)     ? r - 1 : i - 1;
				j = (j == r - 1) ? 0     : j + 1;
			} else {
				++i;
			}
		}

//!# Row exchanges.
		ki = n / 4;
		kj = ki - 1;

		for (i = 0; i < n / 2; i++) {
			if (i != ki) {
				for (j = 0; j < ki; j++) {
					mcswap_var(w, a[(n * i) + j], a[(n * (n / 2 + i)) + j]);
				}
				for (j = 0; j < kj; j++) {
					mcswap_var(w, a[(n * i) + (n - 1 - j)], a[(n * (n / 2 + i)) + (n - 1 - j)]);
				}
			} else {
				for (j = 1; j <= ki; j++) {
					mcswap_var(w, a[(n * ki) + j], a[(n * (n / 2 + ki)) + j]);
				}
				for (j = 0; j < kj; j++) {
					mcswap_var(w, a[(n * ki) + (n - 1 - j)], a[(n * (n / 2 + ki)) + (n - 1 - j)]);
				}
			}
		}
		return 0;
	} else if (n > 2 && ((n % 4) == 0)) {
//!# Doubly-even order.
		r = n / 4;
		for (i = 0; i < r; i++) {
			for (j = 0; j < r; j++) {
				ki = i * 4;
				kj = j * 4;
				for (k = 0; k < 4; k++) {
					for (l = 0; l < 4; l++) {
						a[(n * (ki + k)) + (kj + l)] = mc_cast(float, ((k == l || k + l == 3)
							? (n * n) - (n * (ki + k) + (kj + l))
							: 1       +  n * (ki + k) + (kj + l)
						));
					}
				}
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_magicnxn(const int n, double * a)
{
	int i, j, k, ki, kj;
	int r, l;

	double w;

	if (n > 2 && ((n % 2) == 1)) {
//!# Odd order.
		mc_zerosnxn(n, a);
		i              = 0;
		j              = n / 2;
		k              = 2;
		ki             = i;
		kj             = j;
		a[(n * i) + j] = 1.0;
		for (; k <= n * n; k++) {
			i = i - 1 < 0  ? n - 1 : i - 1;
			j = j + 1 == n ? 0     : j + 1;
			if (a[(n * i) + j] != 0.0) {
				i = ki;
				j = kj;
				i = i + 1 == n ? 0 : i + 1;
			}
			a[(n * i) + j] = mc_cast(double, k);
			ki             = i;
			kj             = j; 
		}
		return 0;
	} else if (n > 2 && ((n % 2 == 0 && (n /2 ) % 2 == 1))) { 
//!# Single-even order.
		r = n / 2;
		i = 0;
		j = r / 2;

//!# Odd order partitions.
		for (k = 1; k <= r * r; k++) {
			a[(n * i) + j]             = mc_cast(double, k            );
			a[(n * (i + r)) + (j + r)] = mc_cast(double, k +     r * r);
			a[(n * i) + (j + r)]       = mc_cast(double, k + 2 * r * r);
			a[(n * (i + r)) + j]       = mc_cast(double, k + 3 * r * r);
			if (!(k % r == 0)) {
				i = (i == 0)     ? r - 1 : i - 1;
				j = (j == r - 1) ? 0     : j + 1;
			} else {
				++i;
			}
		}

//!# Row exchanges.
		ki = n / 4;
		kj = ki - 1;

		for (i = 0; i < n / 2; i++) {
			if (i != ki) {
				for (j = 0; j < ki; j++) {
					mcswap_var(w, a[(n * i) + j], a[(n * (n / 2 + i)) + j]);
				}
				for (j = 0; j < kj; j++) {
					mcswap_var(w, a[(n * i) + (n - 1 - j)], a[(n * (n / 2 + i)) + (n - 1 - j)]);
				}
			} else {
				for (j = 1; j <= ki; j++) {
					mcswap_var(w, a[(n * ki) + j], a[(n * (n / 2 + ki)) + j]);
				}
				for (j = 0; j < kj; j++) {
					mcswap_var(w, a[(n * ki) + (n - 1 - j)], a[(n * (n / 2 + ki)) + (n - 1 - j)]);
				}
			}
		}
		return 0;
	} else if (n > 2 && ((n % 4) == 0)) {
//!# Doubly-even order.
		r = n / 4;
		for (i = 0; i < r; i++) {
			for (j = 0; j < r; j++) {
				ki = i * 4;
				kj = j * 4;
				for (k = 0; k < 4; k++) {
					for (l = 0; l < 4; l++) {
						a[(n * (ki + k)) + (kj + l)] = mc_cast(double, ((k == l || k + l == 3)
							? (n * n) - (n * (ki + k) + (kj + l))
							: 1       +  n * (ki + k) + (kj + l)
						));
					}
				}
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_magicnxnl(const int n, long double * a)
{
	int i, j, k, ki, kj;
	int r, l;

	long double w;

	if (n > 2 && ((n % 2) == 1)) {
//!# Odd order.
		mc_zerosnxnl(n, a);
		i              = 0;
		j              = n / 2;
		k              = 2;
		ki             = i;
		kj             = j;
		a[(n * i) + j] = 1.0L;
		for (; k <= n * n; k++) {
			i = i - 1 < 0  ? n - 1 : i - 1;
			j = j + 1 == n ? 0     : j + 1;
			if (a[(n * i) + j] != 0.0L) {
				i = ki;
				j = kj;
				i = i + 1 == n ? 0 : i + 1;
			}
			a[(n * i) + j] = mc_cast(long double, k);
			ki             = i;
			kj             = j; 
		}
		return 0;
	} else if (n > 2 && ((n % 2 == 0 && (n /2 ) % 2 == 1))) { 
//!# Single-even order.
		r = n / 2;
		i = 0;
		j = r / 2;

//!# Odd order partitions.
		for (k = 1; k <= r * r; k++) {
			a[(n * i) + j]             = mc_cast(long double, k            );
			a[(n * (i + r)) + (j + r)] = mc_cast(long double, k +     r * r);
			a[(n * i) + (j + r)]       = mc_cast(long double, k + 2 * r * r);
			a[(n * (i + r)) + j]       = mc_cast(long double, k + 3 * r * r);
			if (!(k % r == 0)) {
				i = (i == 0)     ? r - 1 : i - 1;
				j = (j == r - 1) ? 0     : j + 1;
			} else {
				++i;
			}
		}

//!# Row exchanges.
		ki = n / 4;
		kj = ki - 1;

		for (i = 0; i < n / 2; i++) {
			if (i != ki) {
				for (j = 0; j < ki; j++) {
					mcswap_var(w, a[(n * i) + j], a[(n * (n / 2 + i)) + j]);
				}
				for (j = 0; j < kj; j++) {
					mcswap_var(w, a[(n * i) + (n - 1 - j)], a[(n * (n / 2 + i)) + (n - 1 - j)]);
				}
			} else {
				for (j = 1; j <= ki; j++) {
					mcswap_var(w, a[(n * ki) + j], a[(n * (n / 2 + ki)) + j]);
				}
				for (j = 0; j < kj; j++) {
					mcswap_var(w, a[(n * ki) + (n - 1 - j)], a[(n * (n / 2 + ki)) + (n - 1 - j)]);
				}
			}
		}
		return 0;
	} else if (n > 2 && ((n % 4) == 0)) {
//!# Doubly-even order.
		r = n / 4;
		for (i = 0; i < r; i++) {
			for (j = 0; j < r; j++) {
				ki = i * 4;
				kj = j * 4;
				for (k = 0; k < 4; k++) {
					for (l = 0; l < 4; l++) {
						a[(n * (ki + k)) + (kj + l)] = mc_cast(long double, ((k == l || k + l == 3)
							? (n * n) - (n * (ki + k) + (kj + l))
							: 1       +  n * (ki + k) + (kj + l)
						));
					}
				}
			}
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_MAGICNXN_H */

/* EOF */