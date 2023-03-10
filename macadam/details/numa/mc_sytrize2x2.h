//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sytrize2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_copy2x2.h>

#ifndef MC_SYTRIZE2X2_H
#define MC_SYTRIZE2X2_H

#pragma mark - mc_sytrize2x2 -

MC_TARGET_FUNC void mc_sytrize2x2f(float b[4], const float a[4], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	if (a != b) {
		mc_copy2x2f(b, a);
	}
	if (f == 1) {
		b[2] = b[1];
	} else {
		b[1] = b[2];
	}
}

MC_TARGET_FUNC void mc_sytrize2x2ff(double b[4], const float a[4], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	mc_copy2x2ff(b, a);
	if (f == 1) {
		b[2] = b[1];
	} else {
		b[1] = b[2];
	}
}

MC_TARGET_FUNC void mc_sytrize2x2(double b[4], const double a[4], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	if (a != b) {
		mc_copy2x2(b, a);
	}
	if (f == 1) {
		b[2] = b[1];
	} else {
		b[1] = b[2];
	}
}

MC_TARGET_FUNC void mc_sytrize2x2l(long double b[4], const long double a[4], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	if (a != b) {
		mc_copy2x2l(b, a);
	}
	if (f == 1) {
		b[2] = b[1];
	} else {
		b[1] = b[2];
	}
}

#endif /* !MC_SYTRIZE2X2_H */

/* EOF */