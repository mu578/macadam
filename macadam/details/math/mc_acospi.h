//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_acospi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_acos.h>

#ifndef MC_ACOSPI_H
#define MC_ACOSPI_H

#pragma mark - mc_acospi -

MC_TARGET_FUNC float mc_acospif(const float x)
{
	const float y = mc_acosf(x);
	return y * MCK_KF(MCK_1_PI);
}

MC_TARGET_FUNC double mc_acospi(const double x)
{
	const double y = mc_acos(x);
	return y * MCK_K(MCK_1_PI);
}

MC_TARGET_FUNC long double mc_acospil(const long double x)
{
	const long double y = mc_acosl(x);
	return y * MCK_KL(MCK_1_PI);
}

#endif /* !MC_ACOSPI_H */

/* EOF */