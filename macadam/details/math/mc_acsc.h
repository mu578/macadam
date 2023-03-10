//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_acsc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_asin.h>

#ifndef MC_ACSC_H
#define MC_ACSC_H

#pragma mark - mc_acsc -

MC_TARGET_FUNC float mc_acscf(const float x)
{
	return mc_asinf(1.0f / x);
}

MC_TARGET_FUNC double mc_acsc(const double x)
{
	return mc_asin(1.0 / x);
}

MC_TARGET_FUNC long double mc_acscl(const long double x)
{
	return mc_asinl(1.0L / x);
}

#endif /* !MC_ACSC_H */

/* EOF */