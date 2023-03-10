//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logdiffexp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log1me.h>

#ifndef MC_LOGDIFFEXP_H
#define MC_LOGDIFFEXP_H

#pragma mark - mc_logdiffexp -

MC_TARGET_FUNC float mc_logdiffexpf(const float x, const float y)
{
	return x + mc_log1mef(x - y);
}

MC_TARGET_FUNC double mc_logdiffexp(const double x, const double y)
{
	return x + mc_log1me(x - y);
}

MC_TARGET_FUNC long double mc_logdiffexpl(const long double x, const long double y)
{
	return x + mc_log1mel(x - y);
}

#endif /* !MC_LOGDIFFEXP_H */

/* EOF */