//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_csch.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_sinh.h>

#ifndef MC_CSCH_H
#define MC_CSCH_H

#pragma mark - mc_csch -

MC_TARGET_FUNC float mc_cschf(const float x)
{
#	if MC_TARGET_EMBEDDED
	const float x0 = mc_expf(-x);
	const float x1 = mc_expf(+x);
	return 2.0f / (x1 - x0);
#	else
	return 1.0f / mc_sinhf(x);
#	endif
}

MC_TARGET_FUNC double mc_csch(const double x)
{
#	if MC_TARGET_EMBEDDED
	const double x0 = mc_exp(-x);
	const double x1 = mc_exp(+x);
	return 2.0 / (x1 - x0);
#	else
	return 1.0 / mc_sinh(x);
#	endif
}

MC_TARGET_FUNC long double mc_cschl(const long double x)
{
#	if MC_TARGET_EMBEDDED
	const long double x0 = mc_expl(-x);
	const long double x1 = mc_expl(+x);
	return 2.0L / (x1 - x0);
#	else
	return 1.0L / mc_sinhl(x);
#	endif
}

#endif /* !MC_CSCH_H */

/* EOF */