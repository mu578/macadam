//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_asec.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_acos.h>

#ifndef MC_ASEC_H
#define MC_ASEC_H

#pragma mark - mc_asec -

MC_TARGET_FUNC float mc_asecf(const float x)
{
	return mc_acosf(1.0f / x);
}

MC_TARGET_FUNC double mc_asec(const double x)
{
	return mc_acos(1.0 / x);
}

MC_TARGET_FUNC long double mc_asecl(const long double x)
{
	return mc_acosl(1.0L / x);
}

#endif /* !MC_ASEC_H */

/* EOF */