//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ffrac.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_modf.h>

#ifndef MC_FFRAC_H
#define MC_FFRAC_H

#pragma mark - mc_ffrac -

MC_TARGET_FUNC float mc_ffracf(const float x)
{
	float y = 0.0f;
	return mc_modff(x, &y);
}

MC_TARGET_FUNC double mc_ffrac(const double x)
{
	double y = 0.0;
	return mc_modf(x, &y);
}

MC_TARGET_FUNC long double mc_ffracl(const long double x)
{
	long double y = 0.0L;
	return mc_modfl(x, &y);
}

#endif /* !MC_FFRAC_H */

/* EOF */