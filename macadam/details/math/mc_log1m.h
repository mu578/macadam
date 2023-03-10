//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_log1m.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log1p.h>

#ifndef MC_LOG1M_H
#define MC_LOG1M_H

#pragma mark - mc_log1m -

MC_TARGET_FUNC float mc_log1mf(const float x)
{
	return mc_log1pf(-x);
}

MC_TARGET_FUNC double mc_log1m(const double x)
{
	return mc_log1p(-x);
}

MC_TARGET_FUNC long double mc_log1ml(const long double x)
{
	return mc_log1pl(-x);
}

#endif /* !MC_LOG1M_H */

/* EOF */