//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zlogspace1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zadd.h>
#include <macadam/details/math/mc_zdiv.h>
#include <macadam/details/math/mc_zmul.h>
#include <macadam/details/math/mc_zpow10.h>
#include <macadam/details/math/mc_zsub.h>

#ifndef MC_ZLOGSPACE1XN
#define MC_ZLOGSPACE1XN

#pragma mark - mc_zlogspace1xn -

MC_TARGET_FUNC int mc_zlogspace1xnf(const int n, float * x_r, float * x_i, const float x1_r, const float x1_i, const float x2_r, const float x2_i)
{
//!# Requires x[n] where 1 < n. Draws a logspace: generates a logarithmically spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	float stepr, stepi;
	double steprd, stepid;
	if (n > 0) {
		if (n < 2) {
			mc_zpow10f(&x_r[0], &x_i[0], x2_r, x2_i);
		} else if (n < 3) {
			mc_zpow10f(&x_r[0], &x_i[0], x1_r, x1_i);
			mc_zpow10f(&x_r[1], &x_i[1], x2_r, x2_i);
		} else {
			mc_zsubf(&stepr, &stepi, x2_r, x2_i, x1_r, x1_i);
			if (n < 0x1000001) {
				mc_zdivf(&stepr, &stepi, stepr, stepi, mc_cast(float, (n - 1)), 0.0f);
			} else {
				steprd = mc_cast(double, stepr);
				stepid = mc_cast(double, stepi);
				mc_zdiv(&steprd, &stepid, steprd, stepid, mc_cast(double, (n - 1)), 0.0);
				stepr  = mc_cast(float, steprd);
				stepi  = mc_cast(float, stepid);
			}
			mc_zpow10f(&x_r[0], &x_i[0], x1_r, x1_i);
			mc_zpow10f(&x_r[(n - 1)], &x_i[(n - 1)], x2_r, x2_i);
			for (; i < (n - 1); i++) {
				mc_zmulf(&x_r[i], &x_i[i], mc_cast(float, i), 0.0f, stepr, stepi);
				mc_zaddf(&x_r[i], &x_i[i], x1_r, x1_i, x_r[i], x_i[i]);
				mc_zpow10f(&x_r[i], &x_i[i], x_r[i], x_i[i]);
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_zlogspace1xnff(const int n, double * x_r, double * x_i, const float x1_r, const float x1_i, const float x2_r, const float x2_i)
{
//!# Requires x[n] where 1 < n. Draws a logspace: generates a logarithmically spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	double stepr, stepi, x1rd, x1id, x2rd, x2id;
	if (n > 0) {
		x1rd = mc_cast(double, x1_r);
		x1id = mc_cast(double, x1_i);
		x2rd = mc_cast(double, x2_r);
		x2id = mc_cast(double, x2_i);
		if (n < 2) {
			mc_zpow10(&x_r[0], &x_i[0], x2rd, x2id);
		} else if (n < 3) {
			mc_zpow10(&x_r[0], &x_i[0], x1rd, x1id);
			mc_zpow10(&x_r[1], &x_i[1], x2rd, x2id);
		} else {
			mc_zsub(&stepr, &stepi, x2rd, x2id, x1rd, x1id);
			mc_zdiv(&stepr, &stepi, stepr, stepi, mc_cast(double, (n - 1)), 0.0);
			mc_zpow10(&x_r[0], &x_i[0], x1rd, x1id);
			mc_zpow10(&x_r[(n - 1)], &x_i[(n - 1)], x2rd, x2id);
			for (; i < (n - 1); i++) {
				mc_zmul(&x_r[i], &x_i[i], mc_cast(double, i), 0.0, stepr, stepi);
				mc_zadd(&x_r[i], &x_i[i], x1rd, x1id, x_r[i], x_i[i]);
				mc_zpow10(&x_r[i], &x_i[i], x_r[i], x_i[i]);
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_zlogspace1xn(const int n, double * x_r, double * x_i, const double x1_r, const double x1_i, const double x2_r, const double x2_i)
{
//!# Requires x[n] where 1 < n. Draws a logspace: generates a logarithmically spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	double stepr, stepi;
	if (n > 0) {
		if (n < 2) {
			mc_zpow10(&x_r[0], &x_i[0], x2_r, x2_i);
		} else if (n < 3) {
			mc_zpow10(&x_r[0], &x_i[0], x1_r, x1_i);
			mc_zpow10(&x_r[1], &x_i[1], x2_r, x2_i);
		} else {
			mc_zsub(&stepr, &stepi, x2_r, x2_i, x1_r, x1_i);
			mc_zdiv(&stepr, &stepi, stepr, stepi, mc_cast(double, (n - 1)), 0.0);
			mc_zpow10(&x_r[0], &x_i[0], x1_r, x1_i);
			mc_zpow10(&x_r[(n - 1)], &x_i[(n - 1)], x2_r, x2_i);
			for (; i < (n - 1); i++) {
				mc_zmul(&x_r[i], &x_i[i], mc_cast(double, i), 0.0, stepr, stepi);
				mc_zadd(&x_r[i], &x_i[i], x1_r, x1_i, x_r[i], x_i[i]);
				mc_zpow10(&x_r[i], &x_i[i], x_r[i], x_i[i]);
			}
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_zlogspace1xnl(const int n, long double * x_r, long double * x_i, const long double x1_r, const long double x1_i, const long double x2_r, const long double x2_i)
{
//!# Requires x[n] where 1 < n. Draws a logspace: generates a logarithmically spaced
//!# vector `x`, i.e n points with spacing between points being (x2-x1)/(n-1).
	int i = 1;
	long double stepr, stepi;
	if (n > 0) {
		if (n < 2) {
			mc_zpow10l(&x_r[0], &x_i[0], x2_r, x2_i);
		} else if (n < 3) {
			mc_zpow10l(&x_r[0], &x_i[0], x1_r, x1_i);
			mc_zpow10l(&x_r[1], &x_i[1], x2_r, x2_i);
		} else {
			mc_zsubl(&stepr, &stepi, x2_r, x2_i, x1_r, x1_i);
			mc_zdivl(&stepr, &stepi, stepr, stepi, mc_cast(long double, (n - 1)), 0.0L);
			mc_zpow10l(&x_r[0], &x_i[0], x1_r, x1_i);
			mc_zpow10l(&x_r[(n - 1)], &x_i[(n - 1)], x2_r, x2_i);
			for (; i < (n - 1); i++) {
				mc_zmull(&x_r[i], &x_i[i], mc_cast(long double, i), 0.0L, stepr, stepi);
				mc_zaddl(&x_r[i], &x_i[i], x1_r, x1_i, x_r[i], x_i[i]);
				mc_zpow10l(&x_r[i], &x_i[i], x_r[i], x_i[i]);
			}
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_ZLOGSPACE1XN */

/* EOF */