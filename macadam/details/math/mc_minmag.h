//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_minmag.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_MINMAG_H
#define MC_MINMAG_H

#pragma mark - mc_minmag -

#	if defined(__clang__)
#	pragma clang diagnostic push
#	pragma clang diagnostic ignored "-Wsign-compare"
#	define mc_minmag(a, b) (((a) < (b)) ? (a) : (b))
#	pragma clang diagnostic pop
#	elif defined(__GNUC__)
#	pragma GCC diagnostic push
#	pragma GCC diagnostic ignored "-Wsign-compare"
#	define mc_minmag(a, b) (((a) < (b)) ? (a) : (b))
#	pragma GCC diagnostic pop
#	elif defined(_MSC_VER)
#	pragma warning(push)
#	pragma warning(disable:4018)
#	define mc_minmag(a, b) (((a) < (b)) ? (a) : (b))
#	pragma warning(pop)
#	else
#	define mc_minmag(a, b) (((a) < (b)) ? (a) : (b))
#	endif

#endif /* !MC_MINMAG_H */

/* EOF */