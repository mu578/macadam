//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ortho3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fisnear.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/numa/mc_copy3x3.h>
#include <macadam/details/numa/mc_dotp3x1.h>
#include <macadam/details/numa/mc_l2norm3x1.h>
#include <macadam/details/numa/mc_zeros3x1.h>
#include <macadam/details/numa/mc_zeros3x3.h>

#ifndef MC_ORTHO3X3_H
#define MC_ORTHO3X3_H

#pragma mark - mc_ortho3x3 -

MC_TARGET_FUNC int mc_ortho3x3f(const float a[9], float tol, float q[9], float * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3] if !null.
//!# A and Q may be the same. Forming a ortho-normalized basis Q using
//!# Modified Gram-Schmidt method + a decimeting column step if norm < tol.
//!# If R is not null upper-triangle is formed.
	const int wantr = mc_nonnullptr(r);

	float bnorm, cnorm, dot;

	if (a != q) {
		mc_copy3x3f(q, a);
	}
	if (wantr) {
		mc_zeros3x3f(r);
	}

	if (tol <= 0.0f) {
		tol = mc_fabsf(tol);
	}
	if (tol == 0.0f) {
		tol = MCLIMITS_EPSILONF;
	}
	bnorm = 0.0f;

	cnorm = mc_l2norm3x1f(3, 0, q);
	if (cnorm != 0.0f) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1f(3, 0, q);
			q[0] = 1.0f;
			if (wantr) {
				r[0] = 0.0f;
			}
		} else {
			bnorm = mc_fmaxf(bnorm, cnorm);
			if (wantr) {
				r[0] = cnorm;
			}
			if (!mc_fisnearf(cnorm, 1.0f, 1)) {
				cnorm = 1.0f / cnorm;
				q[0]  = q[0] * cnorm;
				q[3]  = q[3] * cnorm;
				q[6]  = q[6] * cnorm;
			}
		}
	} else {
		q[0] = 1.0f;
		if (wantr) {
			r[0] = 0.0f;
		}
	}

	dot   = mc_dotp3x1f(3, 0, 1, q, q, 1);
	if (wantr) {
		r[1] = dot;
	}
	q[1]  = q[1] - (dot * q[0]);
	q[4]  = q[4] - (dot * q[3]);
	q[7]  = q[7] - (dot * q[6]);

	cnorm = mc_l2norm3x1f(3, 1, q);
	if (cnorm != 0.0f) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1f(3, 1, q);
			q[4] = 1.0f;
			if (wantr) {
				r[4] = 0.0f;
			}
		} else {
			bnorm = mc_fmaxf(bnorm, cnorm);
			if (wantr) {
				r[4] = cnorm;
			}
			if (!mc_fisnearf(cnorm, 1.0f, 1)) {
				cnorm = 1.0f / cnorm;
				q[1]  = q[1] * cnorm;
				q[4]  = q[4] * cnorm;
				q[7]  = q[7] * cnorm;
			}
		}
	} else {
		q[4] = 1.0f;
		if (wantr) {
			r[4] = 0.0f;
		}
	}

	dot   = mc_dotp3x1f(3, 0, 2, q, q, 1);
	if (wantr) {
		r[2] = dot;
	}
	q[2]  = q[2] - (dot * q[0]);
	q[5]  = q[5] - (dot * q[3]);
	q[8]  = q[8] - (dot * q[6]);

	dot   = mc_dotp3x1f(3, 1, 2, q, q, 1);
	if (wantr) {
		r[5] = dot;
	}
	q[2]  = q[2] - (dot * q[1]);
	q[5]  = q[5] - (dot * q[4]);
	q[8]  = q[8] - (dot * q[7]);

	cnorm = mc_l2norm3x1f(3, 2, q);
	if (cnorm != 0.0f) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1f(3, 2, q);
			q[8] = 1.0f;
			if (wantr) {
				r[8] = 0.0f;
			}
		} else {
			if (wantr) {
				r[8] = cnorm;
			}
			if (!mc_fisnearf(cnorm, 1.0f, 1)) {
				cnorm = 1.0f / cnorm;
				q[2]  = q[2] * cnorm;
				q[5]  = q[5] * cnorm;
				q[8]  = q[8] * cnorm;
			}
		}
	} else {
		q[8] = 1.0f;
		if (wantr) {
			r[8] = 0.0f;
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_ortho3x3ff(const float a[9], float tol, double q[9], double * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3] if !null.
//!# Forming a ortho-normalized basis Q using Modified Gram-Schmidt
//!# method + a decimeting column step if norm < tol. If R is not null
//!# upper-triangle is formed.
	const int wantr = mc_nonnullptr(r);

	double bnorm, cnorm, dot, told;

	mc_copy3x3ff(q, a);

	if (wantr) {
		mc_zeros3x3(r);
	}

	if (tol <= 0.0f) {
		tol = mc_fabsf(tol);
	}
	if (tol == 0.0f) {
		tol = MCLIMITS_EPSILONF;
	}
	told  = mc_cast(double, tol);
	bnorm = 0.0;

	cnorm = mc_l2norm3x1(3, 0, q);
	if (cnorm != 0.0) {
		if (cnorm < told * bnorm) {
			mc_zeros3x1(3, 0, q);
			q[0] = 1.0;
			if (wantr) {
				r[0] = 0.0;
			}
		} else {
			bnorm = mc_fmax(bnorm, cnorm);
			if (wantr) {
				r[0] = cnorm;
			}
			if (!mc_fisnear(cnorm, 1.0, 1)) {
				cnorm = 1.0 / cnorm;
				q[0]  = q[0] * cnorm;
				q[3]  = q[3] * cnorm;
				q[6]  = q[6] * cnorm;
			}
		}
	} else {
		q[0] = 1.0;
		if (wantr) {
			r[0] = 0.0;
		}
	}

	dot   = mc_dotp3x1(3, 0, 1, q, q, 1);
	if (wantr) {
		r[1] = dot;
	}
	q[1]  = q[1] - (dot * q[0]);
	q[4]  = q[4] - (dot * q[3]);
	q[7]  = q[7] - (dot * q[6]);

	cnorm = mc_l2norm3x1(3, 1, q);
	if (cnorm != 0.0) {
		if (cnorm < told * bnorm) {
			mc_zeros3x1(3, 1, q);
			q[4] = 1.0;
			if (wantr) {
				r[4] = 0.0;
			}
		} else {
			bnorm = mc_fmax(bnorm, cnorm);
			if (wantr) {
				r[4] = cnorm;
			}
			if (!mc_fisnear(cnorm, 1.0, 1)) {
				cnorm = 1.0 / cnorm;
				q[1]  = q[1] * cnorm;
				q[4]  = q[4] * cnorm;
				q[7]  = q[7] * cnorm;
			}
		}
	} else {
		q[4] = 1.0;
		if (wantr) {
			r[4] = 0.0;
		}
	}

	dot   = mc_dotp3x1(3, 0, 2, q, q, 1);
	if (wantr) {
		r[2] = dot;
	}
	q[2]  = q[2] - (dot * q[0]);
	q[5]  = q[5] - (dot * q[3]);
	q[8]  = q[8] - (dot * q[6]);

	dot   = mc_dotp3x1(3, 1, 2, q, q, 1);
	if (wantr) {
		r[5] = dot;
	}
	q[2]  = q[2] - (dot * q[1]);
	q[5]  = q[5] - (dot * q[4]);
	q[8]  = q[8] - (dot * q[7]);

	cnorm = mc_l2norm3x1(3, 2, q);
	if (cnorm != 0.0) {
		if (cnorm < told * bnorm) {
			mc_zeros3x1(3, 2, q);
			q[8] = 1.0;
			if (wantr) {
				r[8] = 0.0;
			}
		} else {
			if (wantr) {
				r[8] = cnorm;
			}
			if (!mc_fisnear(cnorm, 1.0, 1)) {
				cnorm = 1.0 / cnorm;
				q[2]  = q[2] * cnorm;
				q[5]  = q[5] * cnorm;
				q[8]  = q[8] * cnorm;
			}
		}
	} else {
		q[8] = 1.0;
		if (wantr) {
			r[8] = 0.0;
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_ortho3x3(const double a[9], double tol, double q[9], double * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3] if !null.
//!# A and Q may be the same. Forming a ortho-normalized basis Q using
//!# Modified Gram-Schmidt method + a decimeting column step if norm < tol.
//!# If R is not null upper-triangle is formed.
	const int wantr = mc_nonnullptr(r);

	double bnorm, cnorm, dot;

	if (a != q) {
		mc_copy3x3(q, a);
	}
	if (wantr) {
		mc_zeros3x3(r);
	}

	if (tol <= 0.0) {
		tol = mc_fabs(tol);
	}
	if (tol == 0.0) {
		tol = MCLIMITS_EPSILON;
	}
	bnorm = 0.0;

	cnorm = mc_l2norm3x1(3, 0, q);
	if (cnorm != 0.0) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1(3, 0, q);
			q[0] = 1.0;
			if (wantr) {
				r[0] = 0.0;
			}
		} else {
			bnorm = mc_fmax(bnorm, cnorm);
			if (wantr) {
				r[0] = cnorm;
			}
			if (!mc_fisnear(cnorm, 1.0, 1)) {
				cnorm = 1.0 / cnorm;
				q[0]  = q[0] * cnorm;
				q[3]  = q[3] * cnorm;
				q[6]  = q[6] * cnorm;
			}
		}
	} else {
		q[0] = 1.0;
		if (wantr) {
			r[0] = 0.0;
		}
	}

	dot   = mc_dotp3x1(3, 0, 1, q, q, 1);
	if (wantr) {
		r[1] = dot;
	}
	q[1]  = q[1] - (dot * q[0]);
	q[4]  = q[4] - (dot * q[3]);
	q[7]  = q[7] - (dot * q[6]);

	cnorm = mc_l2norm3x1(3, 1, q);
	if (cnorm != 0.0) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1(3, 1, q);
			q[4] = 1.0;
			if (wantr) {
				r[4] = 0.0;
			}
		} else {
			bnorm = mc_fmax(bnorm, cnorm);
			if (wantr) {
				r[4] = cnorm;
			}
			if (!mc_fisnear(cnorm, 1.0, 1)) {
				cnorm = 1.0 / cnorm;
				q[1]  = q[1] * cnorm;
				q[4]  = q[4] * cnorm;
				q[7]  = q[7] * cnorm;
			}
		}
	} else {
		q[4] = 1.0;
		if (wantr) {
			r[4] = 0.0;
		}
	}

	dot   = mc_dotp3x1(3, 0, 2, q, q, 1);
	if (wantr) {
		r[2] = dot;
	}
	q[2]  = q[2] - (dot * q[0]);
	q[5]  = q[5] - (dot * q[3]);
	q[8]  = q[8] - (dot * q[6]);

	dot   = mc_dotp3x1(3, 1, 2, q, q, 1);
	if (wantr) {
		r[5] = dot;
	}
	q[2]  = q[2] - (dot * q[1]);
	q[5]  = q[5] - (dot * q[4]);
	q[8]  = q[8] - (dot * q[7]);

	cnorm = mc_l2norm3x1(3, 2, q);
	if (cnorm != 0.0) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1(3, 2, q);
			q[8] = 1.0;
			if (wantr) {
				r[8] = 0.0;
			}
		} else {
			if (wantr) {
				r[8] = cnorm;
			}
			if (!mc_fisnear(cnorm, 1.0, 1)) {
				cnorm = 1.0 / cnorm;
				q[2]  = q[2] * cnorm;
				q[5]  = q[5] * cnorm;
				q[8]  = q[8] * cnorm;
			}
		}
	} else {
		q[8] = 1.0;
		if (wantr) {
			r[8] = 0.0;
		}
	}
	return 0;
}

MC_TARGET_FUNC int mc_ortho3x3l(const long double a[9], long double tol, long double q[9], long double * MC_TARGET_RESTRICT r)
{
//!# Requires a[3 x 3], q[3 x 3] and r[3 x 3] if !null.
//!# A and Q may be the same. Forming a ortho-normalized basis Q using
//!# Modified Gram-Schmidt method + a decimeting column step if norm < tol.
//!# If R is not null upper-triangle is formed.
	const int wantr = mc_nonnullptr(r);

	long double bnorm, cnorm, dot;

	if (a != q) {
		mc_copy3x3l(q, a);
	}
	if (wantr) {
		mc_zeros3x3l(r);
	}

	if (tol <= 0.0L) {
		tol = mc_fabsl(tol);
	}
	if (tol == 0.0L) {
		tol = MCLIMITS_EPSILONF;
	}
	bnorm = 0.0L;

	cnorm = mc_l2norm3x1l(3, 0, q);
	if (cnorm != 0.0L) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1l(3, 0, q);
			q[0] = 1.0L;
			if (wantr) {
				r[0] = 0.0L;
			}
		} else {
			bnorm = mc_fmaxl(bnorm, cnorm);
			if (wantr) {
				r[0] = cnorm;
			}
			if (!mc_fisnearl(cnorm, 1.0L, 1)) {
				cnorm = 1.0L / cnorm;
				q[0]  = q[0] * cnorm;
				q[3]  = q[3] * cnorm;
				q[6]  = q[6] * cnorm;
			}
		}
	} else {
		q[0] = 1.0L;
		if (wantr) {
			r[0] = 0.0L;
		}
	}

	dot   = mc_dotp3x1l(3, 0, 1, q, q, 1);
	if (wantr) {
		r[1] = dot;
	}
	q[1]  = q[1] - (dot * q[0]);
	q[4]  = q[4] - (dot * q[3]);
	q[7]  = q[7] - (dot * q[6]);

	cnorm = mc_l2norm3x1l(3, 1, q);
	if (cnorm != 0.0L) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1l(3, 1, q);
			q[4] = 1.0L;
			if (wantr) {
				r[4] = 0.0L;
			}
		} else {
			bnorm = mc_fmaxl(bnorm, cnorm);
			if (wantr) {
				r[4] = cnorm;
			}
			if (!mc_fisnearl(cnorm, 1.0L, 1)) {
				cnorm = 1.0L / cnorm;
				q[1]  = q[1] * cnorm;
				q[4]  = q[4] * cnorm;
				q[7]  = q[7] * cnorm;
			}
		}
	} else {
		q[4] = 1.0L;
		if (wantr) {
			r[4] = 0.0L;
		}
	}

	dot   = mc_dotp3x1l(3, 0, 2, q, q, 1);
	if (wantr) {
		r[2] = dot;
	}
	q[2]  = q[2] - (dot * q[0]);
	q[5]  = q[5] - (dot * q[3]);
	q[8]  = q[8] - (dot * q[6]);

	dot   = mc_dotp3x1l(3, 1, 2, q, q, 1);
	if (wantr) {
		r[5] = dot;
	}
	q[2]  = q[2] - (dot * q[1]);
	q[5]  = q[5] - (dot * q[4]);
	q[8]  = q[8] - (dot * q[7]);

	cnorm = mc_l2norm3x1l(3, 2, q);
	if (cnorm != 0.0L) {
		if (cnorm < tol * bnorm) {
			mc_zeros3x1l(3, 2, q);
			q[8] = 1.0L;
			if (wantr) {
				r[8] = 0.0L;
			}
		} else {
			if (wantr) {
				r[8] = cnorm;
			}
			if (!mc_fisnearl(cnorm, 1.0L, 1)) {
				cnorm = 1.0L / cnorm;
				q[2]  = q[2] * cnorm;
				q[5]  = q[5] * cnorm;
				q[8]  = q[8] * cnorm;
			}
		}
	} else {
		q[8] = 1.0L;
		if (wantr) {
			r[8] = 0.0L;
		}
	}
	return 0;
}

#endif /* !MC_ORTHO3X3_H */

/* EOF */