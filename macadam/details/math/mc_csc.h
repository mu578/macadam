//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_csc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sin.h>

#ifndef MC_CSC_H
#define MC_CSC_H

#pragma mark - mc_csc -

MC_TARGET_FUNC float mc_cscf(const float x)
{
	return 1.0f / mc_sinf(x);
}

MC_TARGET_FUNC double mc_csc(const double x)
{
	return 1.0 / mc_sin(x);
}

MC_TARGET_FUNC long double mc_cscl(const long double x)
{
	return 1.0L / mc_sinl(x);
}

#endif /* !MC_CSC_H */

/* EOF */