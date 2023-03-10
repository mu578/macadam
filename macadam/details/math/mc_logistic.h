//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logistic.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>

#ifndef MC_LOGISTIC_H
#define MC_LOGISTIC_H

#pragma mark - mc_logistic -

MC_TARGET_FUNC float mc_logisticf(const float x, const float l)
{
	return 1.0f / (mc_expf(-l * x) + 1.0f);
}

MC_TARGET_FUNC double mc_logistic(const double x, const double l)
{
	return 1.0 / (mc_exp(-l * x) + 1.0);
}

MC_TARGET_FUNC long double mc_logisticl(const long double x, const long double l)
{
	return 1.0L / (mc_expl(-l * x) + 1.0L);
}

#endif /* !MC_LOGISTIC_H */

/* EOF */