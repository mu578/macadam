//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_pow2i.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp2i.h>

#ifndef MC_POW2I_H
#define MC_POW2I_H

#pragma mark - mc_pow2i -

MC_TARGET_FUNC float mc_pow2if(const int x)
{
	return mc_exp2if(x);
}

MC_TARGET_FUNC double mc_pow2i(const int x)
{
	return mc_exp2i(x);
}

MC_TARGET_FUNC long double mc_pow2il(const int x)
{
	return mc_exp2il(x);
}

#endif /* !MC_POW2I_H */

/* EOF */