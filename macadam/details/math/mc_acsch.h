//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_acsch.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_asinh.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_ACSCH_H
#define MC_ACSCH_H

#pragma mark - mc_acsch -

MC_TARGET_FUNC float mc_acschf(const float x)
{
#	if MC_TARGET_EMBEDDED
	return mc_logf(1.0f / x + mc_sqrtf(1.0f + mc_raise2f(x)) / mc_fabsf(x));
#	else
	return mc_asinhf(1.0f / x);
#	endif
}

MC_TARGET_FUNC double mc_acsch(const double x)
{
#	if MC_TARGET_EMBEDDED
	return mc_log(1.0 / x + mc_sqrt(1.0 + mc_raise2(x)) / mc_fabs(x));
#	else
	return mc_asinh(1.0 / x);
#	endif
}

MC_TARGET_FUNC long double mc_acschl(const long double x)
{
#	if MC_TARGET_EMBEDDED
	return mc_logl(1.0L / x + mc_sqrtl(1.0L + mc_raise2l(x)) / mc_fabsl(x));
#	else
	return mc_asinhl(1.0L / x);
#	endif
}

#endif /* !MC_ACSCH_H */

/* EOF */