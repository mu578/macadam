//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zmod.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_znorm.h>

#ifndef MC_ZMOD_H
#define MC_ZMOD_H

#pragma mark - mc_zmod -

MC_TARGET_PROC float mc_zmodf(const float z_r, const float z_i)
{
//!# \note: mc_zmod is similar to mc_zabs with additional boundary checks.
	if (mc_isinf(z_r)) {
		return MCK_INFP;
	}
	if (mc_isinf(z_i)) {
		return MCK_INFP;
	}
	const float n = mc_znormf(z_r, z_i);
	return mc_sqrtf(n);
}

MC_TARGET_PROC double mc_zmod(const double z_r, const double z_i)
{
//!# \note: mc_zmod is similar to mc_zabs with additional boundary checks.
	if (mc_isinf(z_r)) {
		return MCK_INFP;
	}
	if (mc_isinf(z_i)) {
		return MCK_INFP;
	}
	const double n = mc_znorm(z_r, z_i);
	return mc_sqrt(n);
}

MC_TARGET_PROC long double mc_zmodl(const long double z_r, const long double z_i)
{
//!# \note: mc_zmod is similar to mc_zabs with additional boundary checks.
	if (mc_isinf(z_r)) {
		return MCK_INFP;
	}
	if (mc_isinf(z_i)) {
		return MCK_INFP;
	}
	const long double n = mc_znorml(z_r, z_i);
	return mc_sqrtl(n);
}

#endif /* !MC_ZMOD_H */

/* EOF */