//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_geometric.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_log1m.h>
#include <macadam/details/rand/mc_randu.h>

#ifndef MC_RAND_GEOMETRIC_H
#define MC_RAND_GEOMETRIC_H

#pragma mark - mc_rand_geometric -

MC_TARGET_FUNC float mc_rand_geometricf(const float p)
{
	const float r =  mc_randuf();
	const float u = r > 0.0f ? mc_logf(r) : 0.5f;
	const float v = 1.0f / mc_log1mf(p);
	return u * v;
}

MC_TARGET_FUNC double mc_rand_geometricff(const float p)
{
	const double r =  mc_randu();
	const double u = r > 0.0 ? mc_log(r) : 0.5;
	const double v = 1.0 / mc_log1m(mc_cast(const double, p));
	return u * v;
}

MC_TARGET_FUNC double mc_rand_geometric(const double p)
{
	const double r =  mc_randu();
	const double u = r > 0.0 ? mc_log(r) : 0.5;
	const double v = 1.0 / mc_log1m(p);
	return u * v;
}

MC_TARGET_FUNC long double mc_rand_geometricl(const long double p)
{
	const long double r = mc_randul();
	const long double u = r > 0.0L ? mc_logl(r) : 0.5L;
	const long double v = 1.0L / mc_log1ml(p);
	return u * v;
}

#endif /* !MC_RAND_GEOMETRIC_H */

/* EOF */