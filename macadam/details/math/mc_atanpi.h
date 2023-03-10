//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_atanpi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_atan.h>

#ifndef MC_ATANPI_H
#define MC_ATANPI_H

#pragma mark - mc_atanpi -

MC_TARGET_FUNC float mc_atanpif(const float x)
{
	const float y = mc_atanf(x);
	return y * MCK_KF(MCK_1_PI);
}

MC_TARGET_FUNC double mc_atanpi(const double x)
{
	const double y = mc_atan(x);
	return y * MCK_K(MCK_1_PI);
}

MC_TARGET_FUNC long double mc_atanpil(const long double x)
{
	const long double y = mc_atanl(x);
	return y * MCK_KL(MCK_1_PI);
}

#endif /* !MC_ATANPI_H */

/* EOF */