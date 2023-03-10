//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lognn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_factorial.h>
#include <macadam/details/math/mc_gammaln.h>
#include <macadam/details/math/mc_log.h>

#ifndef MC_LOGNN_H
#define MC_LOGNN_H

#pragma mark - mc_lognn -

MC_TARGET_FUNC float mc_lognnf(const unsigned int n)
{
//!# Returns log(n!).
	const unsigned int factorial_max = 35U;
	if (n < factorial_max) {
		return mc_logf(mc_factorialf(n));
	}
	return mc_gammalnf(mc_cast(const float, n) + 1.0f);
}

MC_TARGET_FUNC double mc_lognn(const unsigned int n)
{
//!# Returns log(n!).
	const unsigned int factorial_max = 171U;
	if (n < factorial_max) {
		return mc_log(mc_factorial(n));
	}
	return mc_gammaln(mc_cast(const double, n) + 1.0);
}

MC_TARGET_FUNC long double mc_lognnl(const unsigned int n)
{
//!# Returns log(n!).
#	if MC_TARGET_HAVE_LONG_DOUBLE
	const unsigned int factorial_max = 1755U;
#	else
	const unsigned int factorial_max = 171U;
#	endif
	if (n < factorial_max) {
		return mc_logl(mc_factoriall(n));
	}
	return mc_gammalnl(mc_cast(const long double, n) + 1.0L);
}

#endif /* !MC_LOGNN_H */

/* EOF */