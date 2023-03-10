//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_huber_weight.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_huber_loss.h>

#ifndef MC_HUBER_WEIGHT_H
#define MC_HUBER_WEIGHT_H

#pragma mark - mc_huber_weight -

MC_TARGET_FUNC float mc_huber_weightf(const float r)
{
	const float c = 1.345f;
	const float s = 1.0f;
	return mc_huber_lossf(r, c, s, 0);
}

MC_TARGET_FUNC double mc_huber_weight(const double r)
{
	const double c = 1.345;
	const double s = 1.0;
	return mc_huber_loss(r, c, s, 0);
}

MC_TARGET_FUNC long double mc_huber_weightl(const long double r)
{
	const long double c = 1.345L;
	const long double s = 1.0L;
	return mc_huber_lossl(r, c, s, 0);
}

#endif /* !MC_HUBER_WEIGHT_H */

/* EOF */