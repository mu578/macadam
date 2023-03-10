//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_invlogit.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_expit.h>

#ifndef MC_INVLOGIT_H
#define MC_INVLOGIT_H

#pragma mark - mc_invlogit -

MC_TARGET_FUNC float mc_invlogitf(const float x)
{
	return mc_expitf(x);
}

MC_TARGET_FUNC double mc_invlogit(const double x)
{
	return mc_expit(x);
}

MC_TARGET_FUNC long double mc_invlogitl(const long double x)
{
	return mc_expitl(x);
}

#endif /* !MC_INVLOGIT_H */

/* EOF */