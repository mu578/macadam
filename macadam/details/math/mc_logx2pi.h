//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logx2pi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log.h>

#ifndef MC_LOGX2PI_H
#define MC_LOGX2PI_H

#pragma mark - mc_logx2pi -

MC_TARGET_FUNC float mc_logx2pif(const float x)
{
	return MCK_KF(MCK_LOGE2PI) + mc_logf(x);
}

MC_TARGET_FUNC double mc_logx2pi(const double x)
{
	return MCK_K(MCK_LOGE2PI) + mc_log(x);
}

MC_TARGET_FUNC long double mc_logx2pil(const long double x)
{
	return MCK_KL(MCK_LOGE2PI) + mc_logl(x);
}

#endif /* !MC_LOGX2PI_H */

/* EOF */