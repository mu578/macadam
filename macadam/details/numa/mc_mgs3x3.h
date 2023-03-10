//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mgs3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_ortho3x3.h>

#ifndef MC_MGS3X3_H
#define MC_MGS3X3_H

#pragma mark - mc_mgs3x3 -

MC_TARGET_FUNC int mc_mgs3x3f(const float a[9], float q[9], float * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3].
//!# A and Q may be the same. Performing a QR decomposition
//!# of a square matrix using Modified Gram-Schmidt method.
	return mc_ortho3x3f(a, MCLIMITS_TINYF, q, r);
}

MC_TARGET_FUNC int mc_mgs3x3ff(const float a[9], double q[9], double * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3].
//!# Performing a QR decomposition of a square matrix using
//!# Modified Gram-Schmidt method.
	return mc_ortho3x3ff(a, MCLIMITS_TINYF, q, r);
}

MC_TARGET_FUNC int mc_mgs3x3(const double a[9], double q[9], double * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3].
//!# A and Q may be the same. Performing a QR decomposition
//!# of a square matrix using Modified Gram-Schmidt method.
	return mc_ortho3x3(a, MCLIMITS_TINY, q, r);
}

MC_TARGET_FUNC int mc_mgs3x3l(const long double a[9], long double q[9], long double * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3].
//!# A and Q may be the same. Performing a QR decomposition
//!# of a square matrix using Modified Gram-Schmidt method.
	return mc_ortho3x3l(a, MCLIMITS_TINYL, q, r);
}

#endif /* !MC_MGS3X3_H */

/* EOF */