//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logit.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_log1p.h>

#ifndef MC_LOGIT_H
#define MC_LOGIT_H

#pragma mark - mc_logit -

MC_TARGET_FUNC float mc_logitf(const float x)
{
	if (x < MCK_KF(MCK_1_3) || x > MCK_KF(MCK_2_3)) {
		return mc_logf(x / (1.0f - x));
	}
	return mc_log1pf((2.0f * x - 1.0f) / (1.0f - x));
}

MC_TARGET_FUNC double mc_logit(const double x)
{
	if (x < MCK_K(MCK_1_3) || x > MCK_K(MCK_2_3)) {
		return mc_log(x / (1.0 - x));
	}
	return mc_log1p((2.0 * x - 1.0) / (1.0 - x));
}

MC_TARGET_FUNC long double mc_logitl(const long double x)
{
	if (x < MCK_KL(MCK_1_3) || x > MCK_KL(MCK_2_3)) {
		return mc_logl(x / (1.0L - x));
	}
	return mc_log1pl((2.0L * x - 1.0L) / (1.0L - x));
}

#endif /* !MC_LOGIT_H */

/* EOF */