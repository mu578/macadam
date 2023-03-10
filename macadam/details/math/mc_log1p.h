//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_log1p.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_logp1.h>

#ifndef MC_LOG1P_H
#define MC_LOG1P_H

#pragma mark - mc_log1p -

MC_TARGET_FUNC float mc_log1pf(const float x)
{
	return mc_logp1f(x);
}

MC_TARGET_FUNC double mc_log1p(const double x)
{
	return mc_logp1(x);
}

MC_TARGET_FUNC long double mc_log1pl(const long double x)
{
	return mc_logp1l(x);
}

#endif /* !MC_LOG1P_H */

/* EOF */