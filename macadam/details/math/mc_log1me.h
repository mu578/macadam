//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_log1me.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_expm1.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_log1p.h>

#ifndef MC_LOG1ME_H
#define MC_LOG1ME_H

#pragma mark - mc_log1me -

MC_TARGET_FUNC float mc_log1mef(const float x)
{
	return (x > MCK_KF(MCK_LOGE2)) ? mc_log1pf(-mc_expf(-x)) : mc_logf(-mc_expm1f(-x));
}

MC_TARGET_FUNC double mc_log1me(const double x)
{
	return (x > MCK_K(MCK_LOGE2)) ? mc_log1p(-mc_exp(-x)) : mc_log(-mc_expm1(-x));
}

MC_TARGET_FUNC long double mc_log1mel(const long double x)
{
	return (x > MCK_KL(MCK_LOGE2)) ? mc_log1pl(-mc_expl(-x)) : mc_logl(-mc_expm1l(-x));
}

#endif /* !MC_LOG1ME_H */

/* EOF */