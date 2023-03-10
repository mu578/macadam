//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fhrt.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fisneg.h>
#include <macadam/details/math/mc_pow.h>

#ifndef MC_FHRT_H
#define MC_FHRT_H

#pragma mark - mc_fhrt -

MC_TARGET_FUNC float mc_fhrtf(const float x)
{
	if (mc_fisnegf(x)) {
		return -mc_powf(-x, 0.25f);
	}
	return mc_powf(x, 0.25f);
}

MC_TARGET_FUNC double mc_fhrt(const double x)
{
	if (mc_fisneg(x)) {
		return -mc_pow(-x, 0.25);
	}
	return mc_pow(x, 0.25);
}

MC_TARGET_FUNC long double mc_fhrtl(const long double x)
{
	if (mc_fisnegl(x)) {
		return -mc_powl(-x, 0.25L);
	}
	return mc_powl(x, 0.25L);
}

#endif /* !MC_FHRT_H */

/* EOF */