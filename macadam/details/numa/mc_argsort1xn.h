//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_argsort1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcswap.h>

#ifndef MC_ARGSORT1XN_H
#define MC_ARGSORT1XN_H

#pragma mark - mc_argsort1xn_type -

#	define mc_argsort1xn_type(integer_type, n, x, k, f)                                                                                                  \
	mc_scope_begin                                                                                                                                       \
		int __mc_argsort1xn_type_i = 0, __mc_argsort1xn_type_j;                                                                                           \
		integer_type __mc_argsort1xn_type_w;                                                                                                              \
		for (; __mc_argsort1xn_type_i < mc_cast_expr(const int, n); __mc_argsort1xn_type_i++) {                                                           \
			k[__mc_argsort1xn_type_i] = mc_cast(integer_type, __mc_argsort1xn_type_i);                                                                     \
		}                                                                                                                                                 \
		for (__mc_argsort1xn_type_i = 0; __mc_argsort1xn_type_i < mc_cast_expr(const int, n); __mc_argsort1xn_type_i++) {                                 \
			for (__mc_argsort1xn_type_j = 1; __mc_argsort1xn_type_j < (mc_cast_expr(const int, n) - __mc_argsort1xn_type_i); __mc_argsort1xn_type_j++) {   \
				if ((f == 1) ? (x[__mc_argsort1xn_type_j - 1] < x[__mc_argsort1xn_type_j]) : (x[__mc_argsort1xn_type_j - 1] > x[__mc_argsort1xn_type_j])) { \
					mcswap_var(__mc_argsort1xn_type_w, k[__mc_argsort1xn_type_j - 1], k[__mc_argsort1xn_type_j]);                                            \
				}                                                                                                                                           \
			}                                                                                                                                              \
		}                                                                                                                                                 \
	mc_scope_end

#pragma mark - mc_argsort1xn -

MC_TARGET_FUNC void mc_argsort1xnf(const int n, const float * x, int * k, const int f)
{
	mc_argsort1xn_type(int, n, x, k, f);
}

MC_TARGET_FUNC void mc_argsort1xn(const int n, const double * x, int * k, const int f)
{
	mc_argsort1xn_type(int, n, x, k, f);
}

MC_TARGET_FUNC void mc_argsort1xnl(const int n, const long double * x, int * k, const int f)
{
	mc_argsort1xn_type(int, n, x, k, f);
}

#endif /* !MC_ARGSORT1XN_H */

/* EOF */