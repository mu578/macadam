//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cot.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_tan.h>

#ifndef MC_COT_H
#define MC_COT_H

#pragma mark - mc_cot -

MC_TARGET_FUNC float mc_cotf(const float x)
{
	return mc_isfinite(x) && x != 0.0f ? 1.0f / mc_tanf(x) : x;
}

MC_TARGET_FUNC double mc_cot(const double x)
{
	return mc_isfinite(x) && x != 0.0 ? 1.0 / mc_tan(x) : x;
}

MC_TARGET_FUNC long double mc_cotl(const long double x)
{
	return mc_isfinite(x) && x != 0.0L ? 1.0L / mc_tanl(x) : x;
}

#endif /* !MC_COT_H */

/* EOF */