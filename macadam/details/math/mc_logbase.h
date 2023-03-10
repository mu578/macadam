//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_logbase.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_log2.h>
#include <macadam/details/math/mc_log10.h>

#ifndef MC_LOGBASE_H
#define MC_LOGBASE_H

#pragma mark - mc_logbase -

MC_TARGET_FUNC float mc_logbasef(const float x, const int b)
{
//!# Log base function. Computing log10(x) / log10(b=base)
//!# where b > 0. If b is `zero` or `one` returning loge(x).

	float l = MCK_NAN;
	switch (b)
	{
		case 0:
		case 1:
		{
			l = mc_logf(x);
		}
		break;
		case 2:
		{
			l = mc_log2f(x);
		}
		break;
		case 10:
		{
			l = mc_log10f(x);
		}
		break;
		default:
		{
			if (b > 0) {
				l = mc_log10f(x) / mc_log10f(b);
			}
		}
		break;
	}
	return l;
}

MC_TARGET_FUNC double mc_logbase(const double x, const int b)
{
//!# Log base function. Computing log10(x) / log10(b=base)
//!# where b > 0. If b is `zero` or `one` returning loge(x).

	double l = MCK_NAN;
	switch (b)
	{
		case 0:
		case 1:
		{
			l = mc_log(x);
		}
		break;
		case 2:
		{
			l = mc_log2(x);
		}
		break;
		case 10:
		{
			l = mc_log10(x);
		}
		break;
		default:
		{
			if (b > 0) {
				l = mc_log10(x) / mc_log10(b);
			}
		}
		break;
	}
	return l;
}

MC_TARGET_FUNC long double mc_logbasel(const long double x, const int b)
{
//!# Log base function. Computing log10(x) / log10(b=base)
//!# where b > 0. If b is `zero` or `one` returning loge(x).

	long double l = MCK_NAN;
	switch (b)
	{
		case 0:
		case 1:
		{
			l = mc_logl(x);
		}
		break;
		case 2:
		{
			l = mc_log2l(x);
		}
		break;
		case 10:
		{
			l = mc_log10l(x);
		}
		break;
		default:
		{
			if (b > 0) {
				l = mc_log10l(x) / mc_log10l(b);
			}
		}
		break;
	}
	return l;
}

#endif /* !MC_LOGBASE_H */

/* EOF */