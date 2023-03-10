//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fisnint.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fisneg.h>
#include <macadam/details/math/mc_ffrac.h>
#include <macadam/details/math/mc_trunc.h>

#ifndef MC_FISNINT_H
#define MC_FISNINT_H

#pragma mark - mc_fisnint -

MC_TARGET_FUNC int mc_fisnintf(const float x)
{
//!# Determining if `x` is an negative integer value; negative-zero being
//!# a singularity which cannot be represented as an integral type.
	if (x == mc_truncf(x)) {
		return x == 0.0f ? !mc_fisnegf(x) : mc_fisnegf(x);
	}
	if (0.0f == mc_ffracf(x)) {
		return x == 0.0f ? !mc_fisnegf(x) : mc_fisnegf(x);
	}
	return 0;
}

MC_TARGET_FUNC int mc_fisnint(const double x)
{
//!# Determining if `x` is an negative integer value; negative-zero being
//!# a singularity which cannot be represented as an integral type.
	if (x == mc_trunc(x)) {
		return x == 0.0 ? !mc_fisneg(x) : mc_fisneg(x);
	}
	if (0.0 == mc_ffrac(x)) {
		return x == 0.0 ? !mc_fisneg(x) : mc_fisneg(x);
	}
	return 0;
}

MC_TARGET_FUNC int mc_fisnintl(const long double x)
{
//!# Determining if `x` is an negative integer value; negative-zero being
//!# a singularity which cannot be represented as an integral type.
	if (x == mc_truncl(x)) {
		return x == 0.0L ? !mc_fisnegl(x) : mc_fisnegl(x);
	}
	if (0.0L == mc_ffracl(x)) {
		return x == 0.0L ? !mc_fisnegl(x) : mc_fisnegl(x);
	}
	return 0;
}

#endif /* !MC_FISNINT_H */

/* EOF */