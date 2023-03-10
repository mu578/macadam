//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fisint.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fisneg.h>
#include <macadam/details/math/mc_ffrac.h>
#include <macadam/details/math/mc_trunc.h>

#ifndef MC_FISINT_H
#define MC_FISINT_H

#pragma mark - mc_fisint -

MC_TARGET_FUNC int mc_fisintf(const float x)
{
//!# Determining if `x` is an integer value; negative-zero being a
//!# singularity which cannot be represented as an integral type.
	if (x == mc_truncf(x)) {
		return x == 0.0f ? !mc_fisnegf(x) : 1;
	}
	if (0.0f == mc_ffracf(x)) {
		return x == 0.0f ? !mc_fisnegf(x) : 1;
	}
	return 0;
}

MC_TARGET_FUNC int mc_fisint(const double x)
{
//!# Determining if `x` is an integer value; negative-zero being a
//!# singularity which cannot be represented as an integral type.
	if (x == mc_trunc(x)) {
		return x == 0.0 ? !mc_fisneg(x) : 1;
	}
	if (0.0 == mc_ffrac(x)) {
		return x == 0.0 ? !mc_fisneg(x) : 1;
	}
	return 0;
}

MC_TARGET_FUNC int mc_fisintl(const long double x)
{
//!# Determining if `x` is an integer value; negative-zero being a
//!# singularity which cannot be represented as an integral type.
	if (x == mc_truncl(x)) {
		return x == 0.0L ? !mc_fisnegl(x) : 1;
	}
	if (0.0L == mc_ffracl(x)) {
		return x == 0.0L ? !mc_fisnegl(x) : 1;
	}
	return 0;
}

#endif /* !MC_FISINT_H */

/* EOF */