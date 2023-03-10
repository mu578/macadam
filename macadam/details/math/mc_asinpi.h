//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_asinpi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_asin.h>

#ifndef MC_ASINPI_H
#define MC_ASINPI_H

#pragma mark - mc_asinpi -

MC_TARGET_FUNC float mc_asinpif(const float x)
{
	const float y = mc_asinf(x);
	return y * MCK_KF(MCK_1_PI);
}

MC_TARGET_FUNC double mc_asinpi(const double x)
{
	const double y = mc_asin(x);
	return y * MCK_K(MCK_1_PI);
}

MC_TARGET_FUNC long double mc_asinpil(const long double x)
{
	const long double y = mc_asinl(x);
	return y * MCK_KL(MCK_1_PI);
}

#endif /* !MC_ASINPI_H */

/* EOF */