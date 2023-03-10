//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sec.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cos.h>

#ifndef MC_SEC_H
#define MC_SEC_H

#pragma mark - mc_sec -

MC_TARGET_FUNC float mc_secf(const float x)
{
	return 1.0f / mc_cosf(x);
}

MC_TARGET_FUNC double mc_sec(const double x)
{
	return 1.0 / mc_cos(x);
}

MC_TARGET_FUNC long double mc_secl(const long double x)
{
	return 1.0L / mc_cosl(x);
}

#endif /* !MC_SEC_H */

/* EOF */