//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fisodd.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fmod.h>
#include <macadam/details/math/mc_modf.h>

#ifndef MC_FISODD_H
#define MC_FISODD_H

#pragma mark - mc_fisodd -

MC_TARGET_FUNC int mc_fisoddf(const float x, const int frac)
{
	if (x == 0.0f) {
		return 1;
	}
	float y = x;
//!# Checking and extracting integral-part.
	if (frac == 1) {
		mc_modff(y, &y);
	}
//!# Returns if y is odd-integral.
	if (mc_fmodf(y, 1.0f) == 0.0f && mc_fmodf(y, 2.0f) != 0.0f) {
		return 1;
	}
	return 0;
}

MC_TARGET_FUNC int mc_fisodd(const double x, const int frac)
{
	if (x == 0.0) {
		return 1;
	}
	double y = x;
//!# Checking and extracting integral-part.
	if (frac == 1) {
		mc_modf(y, &y);
	}
//!# Returns if y is odd-integral.
	if (mc_fmod(y, 1.0) == 0.0 && mc_fmod(y, 2.0) != 0.0) {
		return 1;
	}
	return 0;
}

MC_TARGET_FUNC int mc_fisoddl(const long double x, const int frac)
{
	if (x == 0.0L) {
		return 1;
	}
	long double y = x;
//!# Checking and extracting integral-part.
	if (frac == 1) {
		mc_modfl(y, &y);
	}
//!# Returns if y is odd-integral.
	if (mc_fmodl(y, 1.0L) == 0.0 && mc_fmodl(y, 2.0L) != 0.0) {
		return 1;
	}
	return 0;
}

#endif /* !MC_FISODD_H */

/* EOF */