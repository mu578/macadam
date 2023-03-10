//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_coth.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_tanh.h>

#ifndef MC_COTH_H
#define MC_COTH_H

#pragma mark - mc_coth -

MC_TARGET_FUNC float mc_cothf(const float x)
{
#	if MC_TARGET_EMBEDDED
	const float x0 = mc_expf(-x);
	const float x1 = mc_expf(+x);
	return (x1 + x0) / (x1 - x0);
#	else
	return 1.0f / mc_tanhf(x);
#	endif
}

MC_TARGET_FUNC double mc_coth(const double x)
{
#	if MC_TARGET_EMBEDDED
	const double x0 = mc_exp(-x);
	const double x1 = mc_exp(+x);
	return (x1 + x0) / (x1 - x0);
#	else
	return 1.0 / mc_tanh(x);
#	endif
}

MC_TARGET_FUNC long double mc_cothl(const long double x)
{
#	if MC_TARGET_EMBEDDED
	const long double x0 = mc_expl(-x);
	const long double x1 = mc_expl(+x);
	return (x1 + x0) / (x1 - x0);
#	else
	return 1.0L / mc_tanhl(x);
#	endif
}

#endif /* !MC_COTH_H */

/* EOF */