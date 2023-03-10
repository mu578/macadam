//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_beta.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_gammasign.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_lbeta.h>

#ifndef MC_BETA_H
#define MC_BETA_H

#pragma mark - mc_beta -

MC_TARGET_FUNC
float mc_betaf(float x, float y)
{
	float r = MCK_NAN;
	if (mc_isnan(x) || mc_isnan(y)) {
		return r;
	}
	if (mc_isinf(x) || mc_isinf(y)) {
		return r;
	}
	if (x == 0.0f && y == 0.0f) {
		return r;
	}
	if (x == 0.0f || y == 0.0f) {
		return MCK_INFP;
	}
	if (x == 1.0f) {
		return 1.0f / y;
	}
	if (y == 1.0f) {
		return 1.0f / x;
	}
	if (x < y) {
		r = x;
		x = y;
		y = r;
	}
	r = x + y;
	if (r <= 0.0f && mc_fisintf(r)) {
		if (mc_ffracf(x) != 0.0f && mc_ffracf(y) != 0.0f) {
			return 0.0f;
		}
	}
	r = mc_lbetaf(x, y);
	if (mc_isnan(r)) {
		return r;
	}
	if (mc_isinf(r)) {
		return r;
	}
	r = mc_expf(r);
	if (x < 0.0f || y < 0.0f) {
		r = r * ((mc_gammasignf(x) * mc_gammasignf(y)) / mc_gammasignf(x + y));
	}
	return r;
}

MC_TARGET_FUNC
double mc_beta(double x, double y)
{
	double r = MCK_NAN;
	if (mc_isnan(x) || mc_isnan(y)) {
		return r;
	}
	if (mc_isinf(x) || mc_isinf(y)) {
		return r;
	}
	if (x == 0.0 && y == 0.0) {
		return r;
	}
	if (x == 0.0 || y == 0.0) {
		return MCK_INFP;
	}
	if (x == 1.0) {
		return 1.0 / y;
	}
	if (y == 1.0) {
		return 1.0 / x;
	}
	if (x < y) {
		r = x;
		x = y;
		y = r;
	}
	r = x + y;
	if (r <= 0.0 && mc_fisint(r)) {
		if (mc_ffrac(x) != 0.0 && mc_ffrac(y) != 0.0) {
			return 0.0;
		}
	}
	r = mc_lbeta(x, y);
	if (mc_isnan(r)) {
		return r;
	}
	if (mc_isinf(r)) {
		return r;
	}
	r = mc_exp(r);
	if (x < 0.0 || y < 0.0) {
		r = r * ((mc_gammasign(x) * mc_gammasign(y)) / mc_gammasign(x + y));
	}
	return r;
}

MC_TARGET_FUNC
long double mc_betal(long double x, long double y)
{
	long double r = MCK_NAN;
	if (mc_isnan(x) || mc_isnan(y)) {
		return r;
	}
	if (mc_isinf(x) || mc_isinf(y)) {
		return r;
	}
	if (x == 0.0L && y == 0.0L) {
		return r;
	}
	if (x == 0.0L || y == 0.0L) {
		return MCK_INFP;
	}
	if (x == 1.0L) {
		return 1.0L / y;
	}
	if (y == 1.0L) {
		return 1.0L / x;
	}
	if (x < y) {
		r = x;
		x = y;
		y = r;
	}
	r = x + y;
if (r <= 0.0L && mc_fisintl(r)) {
		if (mc_ffracl(x) != 0.0L && mc_ffracl(y) != 0.0L) {
			return 0.0L;
		}
	}
	r = mc_lbetal(x, y);
	if (mc_isnan(r)) {
		return r;
	}
	if (mc_isinf(r)) {
		return r;
	}
	r = mc_expl(r);
	if (x < 0.0L || y < 0.0L) {
		r = r * ((mc_gammasignl(x) * mc_gammasignl(y)) / mc_gammasignl(x + y));
	}
	return r;
}

#endif /* !MC_BETA_H */

/* EOF */