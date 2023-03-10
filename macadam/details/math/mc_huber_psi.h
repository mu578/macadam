//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_huber_psi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>

#ifndef MC_HUBER_PSI_H
#define MC_HUBER_PSI_H

#pragma mark - mc_huber_psic -

MC_TARGET_FUNC float mc_huber_psicf(const float r, const float c)
{
	if (r > c) {
		return r;
	}
	if (mc_fabsf(r) <= c) {
		return r;
	}
	return -c;
}

MC_TARGET_FUNC double mc_huber_psic(const double r, const double c)
{
	if (r > c) {
		return r;
	}
	if (mc_fabs(r) <= c) {
		return r;
	}
	return -c;
}

MC_TARGET_FUNC long double mc_huber_psicl(const long double r, const long double c)
{
	if (r > c) {
		return r;
	}
	if (mc_fabsl(r) <= c) {
		return r;
	}
	return -c;
}

#pragma mark - mc_huber_psi -

MC_TARGET_FUNC float mc_huber_psif(const float r)
{
	const float c = 1.345f;
	return mc_huber_psicf(r, c);
}

MC_TARGET_FUNC double mc_huber_psi(const double r)
{
	const double c = 1.345;
	return mc_huber_psic(r, c);
}

MC_TARGET_FUNC long double mc_huber_psil(const long double r)
{
	const long double c = 1.345L;
	return mc_huber_psicl(r, c);
}

#endif /* !MC_HUBER_PSI_H */

/* EOF */