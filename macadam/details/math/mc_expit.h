//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_expit.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>

#ifndef MC_EXPITT_H
#define MC_EXPITT_H

#pragma mark - mc_expit -

MC_TARGET_FUNC float mc_expitf(const float x)
{
	return (1.0f - 2.0f * MCLIMITS_EPSILONF) / (1.0f + mc_expf(-x)) + MCLIMITS_EPSILONF;
}

MC_TARGET_FUNC double mc_expit(const double x)
{
	return (1.0 - 2.0 * MCLIMITS_EPSILON) / (1.0 + mc_exp(-x)) + MCLIMITS_EPSILON;
}

MC_TARGET_FUNC long double mc_expitl(const long double x)
{
	return (1.0L - 2.0L * MCLIMITS_EPSILONL) / (1.0L + mc_expl(-x)) + MCLIMITS_EPSILONL;
}

#endif /* !MC_EXPITT_H */

/* EOF */