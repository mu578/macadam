//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sytrize3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_copy3x3.h>

#ifndef MC_SYTRIZE3X3_H
#define MC_SYTRIZE3X3_H

#pragma mark - mc_sytrize3x3 -

MC_TARGET_FUNC void mc_sytrize3x3f(float b[9], const float a[9], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	if (a != b) {
		mc_copy3x3f(b, a);
	}
	if (f == 1) {
		b[1] = b[3];
		b[2] = b[6];
		b[5] = b[7];
	} else {
		b[3] = b[1];
		b[6] = b[2];
		b[7] = b[5];
	}
}

MC_TARGET_FUNC void mc_sytrize3x3ff(double b[9], const float a[9], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	mc_copy3x3ff(b, a);
	if (f == 1) {
		b[1] = b[3];
		b[2] = b[6];
		b[5] = b[7];
	} else {
		b[3] = b[1];
		b[6] = b[2];
		b[7] = b[5];
	}
}

MC_TARGET_FUNC void mc_sytrize3x3(double b[9], const double a[9], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	if (a != b) {
		mc_copy3x3(b, a);
	}
	if (f == 1) {
		b[1] = b[3];
		b[2] = b[6];
		b[5] = b[7];
	} else {
		b[3] = b[1];
		b[6] = b[2];
		b[7] = b[5];
	}
}

MC_TARGET_FUNC void mc_sytrize3x3l(long double b[9], const long double a[9], const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	if (a != b) {
		mc_copy3x3l(b, a);
	}
	if (f == 1) {
		b[1] = b[3];
		b[2] = b[6];
		b[5] = b[7];
	} else {
		b[3] = b[1];
		b[6] = b[2];
		b[7] = b[5];
	}
}

#endif /* !MC_SYTRIZE3X3_H */

/* EOF */