//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sech.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_exp.h>

#ifndef MC_SECH_H
#define MC_SECH_H

#pragma mark - mc_sech -

MC_TARGET_FUNC float mc_sechf(const float x)
{
#	if MC_TARGET_EMBEDDED
	const float x0 = mc_expf(-x);
	const float x1 = mc_expf(+x);
	return 2.0f / (x1 + x0);
#	else
	return 1.0f / mc_coshf(x);
#	endif
}

MC_TARGET_FUNC double mc_sech(const double x)
{
#	if MC_TARGET_EMBEDDED
	const double x0 = mc_exp(-x);
	const double x1 = mc_exp(+x);
	return 2.0 / (x1 + x0);
#	else
	return 1.0 / mc_cosh(x);
#	endif
}

MC_TARGET_FUNC long double mc_sechl(const long double x)
{
#	if MC_TARGET_EMBEDDED
	const long double x0 = mc_expl(-x);
	const long double x1 = mc_expl(+x);
	return 2.0L / (x1 + x0);
#	else
	return 1.0L / mc_coshl(x);
#	endif
}

#endif /* !MC_SECH_H */

/* EOF */