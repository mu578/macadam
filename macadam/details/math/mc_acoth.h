//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_acoth.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_atanh.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_log.h>

#ifndef MC_ACOTH_H
#define MC_ACOTH_H

#pragma mark - mc_acoth -

MC_TARGET_FUNC float mc_acothf(const float x)
{
//!# For real values x in the domain − inf < x < −1 and 1 < x < inf.
#	if MC_TARGET_EMBEDDED
		if (mc_fabsf(x) <= 1.0f) {
			return MCK_NAN;
		}
		return 0.5f * mc_logf((x + 1.0f) / (x - 1.0f));
#	else
	return mc_atanhf(1.0f / x);
#	endif
}

MC_TARGET_FUNC double mc_acoth(const double x)
{
//!# For real values x in the domain − inf < x < −1 and 1 < x < inf.
#	if MC_TARGET_EMBEDDED
		if (mc_fabs(x) <= 1.0) {
			return MCK_NAN;
		}
		return 0.5 * mc_log((x + 1.0) / (x - 1.0));
#	else
	return mc_atanh(1.0 / x);
#	endif
}

MC_TARGET_FUNC long double mc_acothl(const long double x)
{
//!# For real values x in the domain − inf < x < −1 and 1 < x < inf.
#	if MC_TARGET_EMBEDDED
		if (mc_fabsl(x) <= 1.0L) {
			return MCK_NAN;
		}
		return 0.5L * mc_logl((x + 1.0L) / (x - 1.0L));
#	else
	return mc_atanhl(1.0L / x);
#	endif
}

#endif /* !MC_ACOTH_H */

/* EOF */