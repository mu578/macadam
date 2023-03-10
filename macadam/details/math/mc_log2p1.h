//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_log2p1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log1p.h>

#ifndef MC_LOG2P1_H
#define MC_LOG2P1_H

#pragma mark - mc_log2p1 -

MC_TARGET_FUNC float mc_log2p1f(const float x)
{
	if (x == 0.0f) {
		return x;
	}
	const float y = x * MCK_KF(MCK_1_LOGE2);
	return mc_log1pf(y);
}

MC_TARGET_FUNC double mc_log2p1(const double x)
{
	if (x == 0.0) {
		return x;
	}
	const double y = x * MCK_K(MCK_1_LOGE2);
	return mc_log1p(y);
}

MC_TARGET_FUNC long double mc_log2p1l(const long double x)
{
	if (x == 0.0L) {
		return x;
	}
#	if (MC_TARGET_C99 || MC_TARGET_CPP17) && defined(M_LN2l)
	const long double y = x / M_LN2l;
#	else
	const long double y = x * MCK_KL(MCK_1_LOGE2);
#	endif
	return mc_log1pl(y);
}

#endif /* !MC_LOG2P1_H */

/* EOF */