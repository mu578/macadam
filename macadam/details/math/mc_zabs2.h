//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zabs2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>

#ifndef MC_ZABS2_H
#define MC_ZABS2_H

#pragma mark - mc_zabs2 -

MC_TARGET_PROC float mc_zabs2f(const float z_r, const float z_i)
{
	return mc_raise2f(z_i) + mc_raise2f(z_r);
}

MC_TARGET_PROC double mc_zabs2(const double z_r, const double z_i)
{
	return mc_raise2(z_i) + mc_raise2(z_r);
}

MC_TARGET_PROC long double mc_zabs2l(const long double z_r, const long double z_i)
{
	return mc_raise2l(z_i) + mc_raise2l(z_r);
}

#endif /* !MC_ZABS2_H */

/* EOF */