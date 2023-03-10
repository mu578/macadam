//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_access.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>
#include <macadam/mcswap.h>

#ifndef MC_BLAS_ACCESS_H
#define MC_BLAS_ACCESS_H

#	define mc_nonblas_vector_at(g, gi) g[(gi)]
#	define mc_nonblas_matrix_rmj_at(g, mg, ng, gi, gj) g[((gi) * (ng)) + (gj)]
#	define mc_nonblas_matrix_cmj_at(g, mg, ng, gi, gj) g[((gj) * (mg)) + (gi)]

#	define mc_blas_vector_at(g, gi) mc_nonblas_vector_at(g, ((gi) - 1))

#	if MC_TARGET_BLAS_USE_CLAYOUT
#		define mc_blas_matrix_at(g, mg, ng, gj, gi) \
			mc_nonblas_matrix_rmj_at(g, mg, ng, ((gj) - 1), ((gi) - 1))
#	else
#		define mc_blas_matrix_at(g, mg, ng, gj, gi) \
			mc_nonblas_matrix_cmj_at(g, mg, ng, ((gj) - 1), ((gi) - 1))
#	endif

#endif /* !MC_BLAS_ACCESS_H */

/* EOF */