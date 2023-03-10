//
// # -*- coding: utf-8, tab-width: 3 -*-

// mcswap.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MCSWAP_H
#define MCSWAP_H

#pragma mark - mcswap_var -

#	define mcswap_var(var, a, b) \
		(var) = (a);              \
		(a)   = (b);              \
		(b)   = (var)

#pragma mark - mcswap_type -

#	define mcswap_type(type, a, b)               \
	mc_scope_begin                               \
		type __mcswap_type_aa = (a);              \
		type __mcswap_type_bb = (b);              \
		(a)                   = __mcswap_type_bb; \
		(b)                   = __mcswap_type_aa; \
	mc_scope_end

#pragma mark - mcswap -

#	if MC_TARGET_CPP98
#	define mcswap(a, b) ::std::swap(a, b)
#	elif MC_TARGET_HAVE_AUTOTYPE
#	define mcswap(a, b)                              \
	__extension__ ({                                 \
		MC_TARGET_AUTOTYPE __mcswap_aa = (a);         \
		MC_TARGET_AUTOTYPE __mcswap_bb = (b);         \
		(a)                            = __mcswap_bb; \
		(b)                            = __mcswap_aa; \
		MC_NULLPTR;                                   \
	})
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcswap(a, b) mcswap_type(MC_TARGET_TYPEOF(a), a, b)
#	else
#	define mcswap(a, b) mc_unused(a); mc_unused(b)
#	endif

#endif /* !MCSWAP_H */

/* EOF */