//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_znorm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_ZNORM_H
#define MC_ZNORM_H

#pragma mark - mc_znorm -

MC_TARGET_PROC float mc_znormf(const float z_r, const float z_i)
{
	if (mc_isinf(z_r)) {
		return MCK_INFP;
	}
	if (mc_isinf(z_i)) {
		return MCK_INFP;
	}
	return mc_raise2f(z_r) + mc_raise2f(z_i);
}

MC_TARGET_PROC double mc_znorm(const double z_r, const double z_i)
{
	if (mc_isinf(z_r)) {
		return MCK_INFP;
	}
	if (mc_isinf(z_i)) {
		return MCK_INFP;
	}
	return mc_raise2(z_r) + mc_raise2(z_i);
}

MC_TARGET_PROC long double mc_znorml(const long double z_r, const long double z_i)
{
	if (mc_isinf(z_r)) {
		return MCK_INFP;
	}
	if (mc_isinf(z_i)) {
		return MCK_INFP;
	}
	return mc_raise2l(z_r) + mc_raise2l(z_i);
}

#endif /* !MC_ZNORM_H */

/* EOF */