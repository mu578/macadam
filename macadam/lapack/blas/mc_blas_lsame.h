//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_lsame.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    lsame returns `one` if `ca` is the same character as `cb` regardless of case otherwise `zero`.

 *
 * \synopsis
 *    int lsame(ca, cb)
 *    char ca, cb
 *
 * \purpose
 *    lsame returns `one` if `ca` is the same character as `cb` regardless of case otherwise `zero`.
 *
 * \parameters
 *    [in] ca - char. The left single characters to be compared.
 *    [in] cb - char. The right single characters to be compared.
 *
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#include <macadam/details/mc_target.h>

#ifndef MC_BLAS_LSAME_H
#define MC_BLAS_LSAME_H

MC_TARGET_FUNC int mc_blas_lsame(const char ca, const char cb)
{
#	if MC_TARGET_CPP98
	return (ca == cb) ? 1 : ::toupper(ca) == ::toupper(cb) ? 1 : 0;
#	else
	return (ca == cb) ? 1 : toupper(ca) == toupper(cb) ? 1 : 0;
#	endif
}

#endif /* !MC_BLAS_LSAME_H */

/* EOF */