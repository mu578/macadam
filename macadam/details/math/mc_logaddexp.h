//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logaddexp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_log1p.h>

#ifndef MC_LOGADDEXP_H
#define MC_LOGADDEXP_H

#pragma mark - mc_logaddexp -

MC_TARGET_FUNC float mc_logaddexpf(const float x, const float y)
{
	return ((x == y)
		? (x + MCK_KF(MCK_LOGE2))
		: ((x - y) > 0 ? x + mc_log1pf(mc_expf(-(x - y))) : y + mc_log1pf(mc_expf((x - y))))
	);
}

MC_TARGET_FUNC double mc_logaddexp(const double x, const double y)
{
	return ((x == y)
		? (x + MCK_K(MCK_LOGE2))
		: ((x - y) > 0 ? x + mc_log1p(mc_exp(-(x - y))) : y + mc_log1p(mc_exp((x - y))))
	);
}

MC_TARGET_FUNC long double mc_logaddexpl(const long double x, const long double y)
{
#	if (MC_TARGET_C99 || MC_TARGET_CPP17) && defined(M_LN2l)
		return ((x == y)
			? (x + M_LN2l)
			: ((x - y) > 0 ? x + mc_log1pl(mc_expl(-(x - y))) : y + mc_log1pl(mc_expl((x - y))))
		);
#	else
		return ((x == y)
			? (x + MCK_K(MCK_LOGE2))
			: ((x - y) > 0 ? x + mc_log1pl(mc_expl(-(x - y))) : y + mc_log1pl(mc_expl((x - y))))
		);
#	endif
}

#endif /* !MC_LOGADDEXP_H */

/* EOF */