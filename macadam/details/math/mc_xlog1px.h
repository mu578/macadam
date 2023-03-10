//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_xlog1px.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_xlog1py.h>

#ifndef MC_XLOG1PX_H
#define MC_XLOG1PX_H

#pragma mark - mc_xlog1px -

MC_TARGET_FUNC float mc_xlog1pxf(const float x)
{
	return mc_xlog1pyf(x, x);
}

MC_TARGET_FUNC double mc_xlog1px(const double x)
{
	return mc_xlog1py(x, x);
}

MC_TARGET_FUNC long double mc_xlog1pxl(const long double x)
{
	return mc_xlog1pyl(x, x);
}

#endif /* !MC_XLOG1PX_H */

/* EOF */