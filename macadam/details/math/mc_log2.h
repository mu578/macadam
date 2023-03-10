//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_log2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_LOG2_H
#define MC_LOG2_H

#pragma mark - mc_log2 -

MC_TARGET_FUNC float mc_log2f(const float x)
{
#	if MC_TARGET_HAVE_LOG2
#	if MC_TARGET_CPP98
	return ::log2f(x);
#	else
	return log2f(x);
#	endif
#	else
#	if MC_TARGET_CPP98
	return ::logf(x) * MCK_KF(MCK_1_LOGE2);
#	else
	return logf(x) * MCK_KF(MCK_1_LOGE2);
#	endif
#	endif
}

MC_TARGET_FUNC double mc_log2(const double x)
{
#	if MC_TARGET_HAVE_LOG2
#	if MC_TARGET_CPP98
	return ::log2(x);
#	else
	return log2(x);
#	endif
#	else
#	if MC_TARGET_CPP98
	return ::log(x) * MCK_K(MCK_1_LOGE2);
#	else
	return log(x) * MCK_K(MCK_1_LOGE2);
#	endif
#	endif
}

MC_TARGET_FUNC long double mc_log2l(const long double x)
{
#	if MC_TARGET_HAVE_LOG2
#	if MC_TARGET_CPP98
	return ::log2l(x);
#	else
	return log2l(x);
#	endif
#	else
#	if MC_TARGET_CPP98
	return ::logl(x) * MCK_KL(MCK_1_LOGE2);
#	else
	return logl(x) * MCK_KL(MCK_1_LOGE2);
#	endif
#	endif
}

#endif /* !MC_LOG1M_H */

/* EOF */