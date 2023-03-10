//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_asech.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_acosh.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_ASECH_H
#define MC_ASECH_H

#pragma mark - mc_asech -

MC_TARGET_FUNC float mc_asechf(const float x)
{
#	if MC_TARGET_EMBEDDED
	if (x <= 0.0f || x > 1.0f) {
		return MCK_NAN;
	}
	return mc_logf((1.0f + mc_sqrtf(1.0f - mc_raise2f(x))) / x);
#	else
	return mc_acoshf(1.0f / x);
#	endif
}

MC_TARGET_FUNC double mc_asech(const double x)
{
#	if MC_TARGET_EMBEDDED
	if (x <= 0.0 || x > 1.0) {
		return MCK_NAN;
	}
	return mc_log((1.0 + mc_sqrt(1.0 - mc_raise2(x))) / x);
#	else
	return mc_acosh(1.0 / x);
#	endif
}

MC_TARGET_FUNC long double mc_asechl(const long double x)
{
#	if MC_TARGET_EMBEDDED
	if (x <= 0.0L || x > 1.0L) {
		return MCK_NAN;
	}
	return mc_logl((1.0L + mc_sqrtl(1.0L - mc_raise2l(x))) / x);
#	else
	return mc_acoshl(1.0L / x);
#	endif
}

#endif /* !MC_ASECH_H */

/* EOF */