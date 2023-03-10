//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fisneg.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>

#ifndef MC_FISNEG_H
#define MC_FISNEG_H

#pragma mark - mc_fisneg -

MC_TARGET_FUNC int mc_fisnegf(const float x)
{
//!# Determining if `x` is a negative real value; negative-zero
//!# being a singularity an additional sign check is performed.
	return x == 0.0f ? (mc_copysignf(1.0f, x) < 0.0f ? 1 : 0) : (x < 0.0f ? 1 : 0);
}

MC_TARGET_FUNC int mc_fisneg(const double x)
{
//!# Determining if `x` is a negative real value; negative-zero
//!# being a singularity an additional sign check is performed.
	return x == 0.0 ? (mc_copysign(1.0, x) < 0.0 ? 1 : 0) : (x < 0.0 ? 1 : 0);
}

MC_TARGET_FUNC int mc_fisnegl(const long double x)
{
//!# Determining if `x` is a negative real value; negative-zero
//!# being a singularity an additional sign check is performed.
	return x == 0.0L ? (mc_copysignl(1.0L, x) < 0.0L ? 1 : 0) : (x < 0.0L ? 1 : 0);
}

#endif /* !MC_FISNEG_H */

/* EOF */