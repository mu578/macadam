//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulabnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_mulab2x2.h>
#include <macadam/details/numa/mc_mulab3x3.h>
#include <macadam/details/numa/mc_mulabmxn.h>

#ifndef MC_MULABNXN_H
#define MC_MULABNXN_H

#pragma mark - mc_mulabnxn -

MC_TARGET_FUNC void mc_mulabnxnf(const int n, float * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# c=a*b
	if (n == 1) {
		c[0] = a[0] * b[0];
	} else if (n == 2) {
		mc_mulab2x2f(c, a, b);
	} else if (n == 3) {
		mc_mulab3x3f(c, a, b);
	} else {
		mc_mulabmxnf(n, n, n, c, a, b);
	}
}

MC_TARGET_FUNC void mc_mulabnxnff(const int n, double * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# c=a*b
	if (n == 1) {
		c[0] = a[0] * b[0];
	} else if (n == 2) {
		mc_mulab2x2ff(c, a, b);
	} else if (n == 3) {
		mc_mulab3x3ff(c, a, b);
	} else {
		mc_mulabmxnff(n, n, n, c, a, b);
	}
}

MC_TARGET_FUNC void mc_mulabnxnfd(const int n, double * MC_TARGET_RESTRICT c, const float * a, const double * b)
{
//!# c=a*b
	if (n == 1) {
		c[0] = a[0] * b[0];
	} else if (n == 2) {
		mc_mulab2x2fd(c, a, b);
	} else if (n == 3) {
		mc_mulab3x3fd(c, a, b);
	} else {
		mc_mulabmxnfd(n, n, n, c, a, b);
	}
}

MC_TARGET_FUNC void mc_mulabnxndf(const int n, double * MC_TARGET_RESTRICT c, const double * a, const float * b)
{
//!# c=a*b
	if (n == 1) {
		c[0] = a[0] * b[0];
	} else if (n == 2) {
		mc_mulab2x2df(c, a, b);
	} else if (n == 3) {
		mc_mulab3x3df(c, a, b);
	} else {
		mc_mulabmxndf(n, n, n, c, a, b);
	}
}

MC_TARGET_FUNC void mc_mulabnxn(const int n, double * MC_TARGET_RESTRICT c, const double * a, const double * b)
{
//!# c=a*b
	if (n == 1) {
		c[0] = a[0] * b[0];
	} else if (n == 2) {
		mc_mulab2x2(c, a, b);
	} else if (n == 3) {
		mc_mulab3x3(c, a, b);
	} else {
		mc_mulabmxn(n, n, n, c, a, b);
	}
}

MC_TARGET_FUNC void mc_mulabnxnl(const int n, long double * MC_TARGET_RESTRICT c, const long double * a, const long double * b)
{
//!# c=a*b
	if (n == 1) {
		c[0] = a[0] * b[0];
	} else if (n == 2) {
		mc_mulab2x2l(c, a, b);
	} else if (n == 3) {
		mc_mulab3x3l(c, a, b);
	} else {
		mc_mulabmxnl(n, n, n, c, a, b);
	}
}

#endif /* !MC_MULABNXN_H */

/* EOF */