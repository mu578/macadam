//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sort1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_copy1xn.h>
#include <macadam/mcswap.h>

#ifndef MC_SORT1XN_H
#define MC_SORT1XN_H

#pragma mark - mc_sort1xn_type -

#	define mc_sort1xn_type(scalar_type, n, x, f)                                                                                                \
	mc_scope_begin                                                                                                                              \
		size_t __mc_sort1xn_type_i = 0, __mc_sort1xn_type_j;                                                                                     \
		scalar_type __mc_sort1xn_type_w;                                                                                                         \
		if ((n) > 0) {                                                                                                                           \
			for (; __mc_sort1xn_type_i < mc_cast_expr(size_t, n); __mc_sort1xn_type_i++) {                                                        \
				for (__mc_sort1xn_type_j = 1; __mc_sort1xn_type_j < (mc_cast_expr(size_t, n) - __mc_sort1xn_type_i); __mc_sort1xn_type_j++) {      \
					if ((f == 1) ? (x[__mc_sort1xn_type_j - 1] < x[__mc_sort1xn_type_j]) : (x[__mc_sort1xn_type_j - 1] > x[__mc_sort1xn_type_j])) { \
						mcswap_var(__mc_sort1xn_type_w, x[__mc_sort1xn_type_j - 1], x[__mc_sort1xn_type_j]);                                         \
					}                                                                                                                               \
				}                                                                                                                                  \
			}                                                                                                                                     \
		}                                                                                                                                        \
	mc_scope_end

#pragma mark - mc_sort1xn -

MC_TARGET_FUNC void mc_sort1xnf(const int n, float * y, const float * x, const int f)
{
	if (x != y) {
		mc_copy1xnf(n, y, x);
	}
	mc_sort1xn_type(float, n, y, f);
}

MC_TARGET_FUNC void mc_sort1xnff(const int n, double * y, const float * x, const int f)
{
	mc_copy1xnff(n, y, x);
	mc_sort1xn_type(double, n, y, f);
}

MC_TARGET_FUNC void mc_sort1xn(const int n, double * y, const double * x, const int f)
{
	if (x != y) {
		mc_copy1xn(n, y, x);
	}
	mc_sort1xn_type(double, n, y, f);
}

MC_TARGET_FUNC void mc_sort1xnl(const int n, long double * y, const long double * x, const int f)
{
	if (x != y) {
		mc_copy1xnl(n, y, x);
	}
	mc_sort1xn_type(long double, n, y, f);
}

#endif /* !MC_SORT1XN_H */

/* EOF */