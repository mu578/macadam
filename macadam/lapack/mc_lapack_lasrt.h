//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasrt.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_LAPACKE_LASRT_H
#define MC_LAPACKE_LASRT_H

#pragma mark - mc_lapack_slasrt -

MC_TARGET_FUNC void mc_lapack_slasrt(const char id, const int n, float * d, int * info)
{
	int select = 20;
	int dir, endd, i, j, start, stkpnt;
	float d1, d2, d3, dmnmx, tmp;

	int stack[2 * 32] = { 0 };

	*info = 0;
	 dir  = -1;
	if (mc_blas_lsame(id, 'D')) {
		dir = 0;
	} else if (mc_blas_lsame(id, 'I')) {
		dir = 1;
	}
	if (dir == -1) {
		*info = -1;
	} else if (n < 0) {
		*info = -2;
	}
	if (*info != 0) {
		mc_blas_xerbla("SLASRT", -(*info));
		return;
	}

	if (n <= 1) {
		return;
	}

	stkpnt                                = 1;
	mc_blas_matrix_at(stack, 2, 32, 1, 1) = 1;
	mc_blas_matrix_at(stack, 2, 32, 2, 1) = n;

F10:
	start  = mc_blas_matrix_at(stack, 2, 32, 1, stkpnt);
	endd   = mc_blas_matrix_at(stack, 2, 32, 2, stkpnt);
	stkpnt = stkpnt - 1;
	if (endd - start <= select && endd - start > 0) {
		if (dir == 0) {
			for (i = (start + 1); i <= endd; ++i) {
				for (j = i; j >= (start + 1); --j) {
					if (mc_blas_vector_at(d, j) > mc_blas_vector_at(d, j - 1)) {
						dmnmx                       = mc_blas_vector_at(d, j);
						mc_blas_vector_at(d, j)     = mc_blas_vector_at(d, j - 1);
						mc_blas_vector_at(d, j - 1) = dmnmx;
					} else {
						goto F30;
					}
				}
F30:;
			}
		} else {
			for (i = (start + 1); i <= endd; ++i) {
				for (j = i; j >= (start + 1); --j) {
					if (mc_blas_vector_at(d, j) < mc_blas_vector_at(d, j - 1)) {
						dmnmx                       = mc_blas_vector_at(d, j);
						mc_blas_vector_at(d, j)     = mc_blas_vector_at(d, j - 1);
						mc_blas_vector_at(d, j - 1) = dmnmx;
					} else {
						goto F50;
					}
				}
F50:;
			}
		}

	} else if (endd - start > select) {
		d1 = mc_blas_vector_at(d, start);
		d2 = mc_blas_vector_at(d, endd);
		i  = (start + endd) / 2;
		d3 = mc_blas_vector_at(d, i);
		if (d1 < d2) {
			if (d3 < d1) {
				dmnmx = d1;
			} else if (d3 < d2) {
				dmnmx = d3;
			} else {
				dmnmx = d2;
			}
		} else {
			if (d3 < d2) {
				dmnmx = d2;
			} else if (d3 < d1) {
				dmnmx = d3;
			} else {
				dmnmx = d1;
			}
		}
		if (dir == 0) {
			i = start - 1;
			j = endd + 1;
F60:
F70:
			j = j - 1;
			if (mc_blas_vector_at(d, j) < dmnmx) {
				goto F70;
			}
F80:
			i = i + 1;
			if (mc_blas_vector_at(d, i) > dmnmx) {
				goto F80;
			}
			if (i < j) {
				tmp                     = mc_blas_vector_at(d, i);
				mc_blas_vector_at(d, i) = mc_blas_vector_at(d, j);
				mc_blas_vector_at(d, j) = tmp;
				goto F60;
			}
			if (j - start > endd - j - 1) {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
			} else {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
			}
		} else {
			i = start - 1;
			j = endd + 1;
F90:
F100:
			j = j - 1;
			if (mc_blas_vector_at(d, j) > dmnmx) {
				goto F100;
			}
	F110:
			i = i + 1;
			if (mc_blas_vector_at(d, i) < dmnmx) {
				goto F110;
			}
			if (i < j) {
				tmp                     = mc_blas_vector_at(d, i);
				mc_blas_vector_at(d, i) = mc_blas_vector_at(d, j);
				mc_blas_vector_at(d, j) = tmp;
				goto F90;
			}
			if (j - start > endd - j - 1) {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
			} else {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
			}
		}
	}
	if (stkpnt > 0) {
		goto F10;
	}
}

#pragma mark - mc_lapack_dlasrt -

MC_TARGET_FUNC void mc_lapack_dlasrt(const char id, const int n, double * d, int * info)
{
	int select = 20;
	int dir, endd, i, j, start, stkpnt;
	double d1, d2, d3, dmnmx, tmp;

	int stack[2 * 32] = { 0 };

	*info = 0;
	 dir  = -1;
	if (mc_blas_lsame(id, 'D')) {
		dir = 0;
	} else if (mc_blas_lsame(id, 'I')) {
		dir = 1;
	}
	if (dir == -1) {
		*info = -1;
	} else if (n < 0) {
		*info = -2;
	}
	if (*info != 0) {
		mc_blas_xerbla("DLASRT", -(*info));
		return;
	}

	if (n <= 1) {
		return;
	}

	stkpnt                                = 1;
	mc_blas_matrix_at(stack, 2, 32, 1, 1) = 1;
	mc_blas_matrix_at(stack, 2, 32, 2, 1) = n;

F10:
	start  = mc_blas_matrix_at(stack, 2, 32, 1, stkpnt);
	endd   = mc_blas_matrix_at(stack, 2, 32, 2, stkpnt);
	stkpnt = stkpnt - 1;
	if (endd - start <= select && endd - start > 0) {
		if (dir == 0) {
			for (i = (start + 1); i <= endd; ++i) {
				for (j = i; j >= (start + 1); --j) {
					if (mc_blas_vector_at(d, j) > mc_blas_vector_at(d, j - 1)) {
						dmnmx                       = mc_blas_vector_at(d, j);
						mc_blas_vector_at(d, j)     = mc_blas_vector_at(d, j - 1);
						mc_blas_vector_at(d, j - 1) = dmnmx;
					} else {
						goto F30;
					}
				}
F30:;
			}
		} else {
			for (i = (start + 1); i <= endd; ++i) {
				for (j = i; j >= (start + 1); --j) {
					if (mc_blas_vector_at(d, j) < mc_blas_vector_at(d, j - 1)) {
						dmnmx                       = mc_blas_vector_at(d, j);
						mc_blas_vector_at(d, j)     = mc_blas_vector_at(d, j - 1);
						mc_blas_vector_at(d, j - 1) = dmnmx;
					} else {
						goto F50;
					}
				}
F50:;
			}
		}

	} else if (endd - start > select) {
		d1 = mc_blas_vector_at(d, start);
		d2 = mc_blas_vector_at(d, endd);
		i  = (start + endd) / 2;
		d3 = mc_blas_vector_at(d, i);
		if (d1 < d2) {
			if (d3 < d1) {
				dmnmx = d1;
			} else if (d3 < d2) {
				dmnmx = d3;
			} else {
				dmnmx = d2;
			}
		} else {
			if (d3 < d2) {
				dmnmx = d2;
			} else if (d3 < d1) {
				dmnmx = d3;
			} else {
				dmnmx = d1;
			}
		}
		if (dir == 0) {
			i = start - 1;
			j = endd + 1;
F60:
F70:
			j = j - 1;
			if (mc_blas_vector_at(d, j) < dmnmx) {
				goto F70;
			}
F80:
			i = i + 1;
			if (mc_blas_vector_at(d, i) > dmnmx) {
				goto F80;
			}
			if (i < j) {
				tmp                     = mc_blas_vector_at(d, i);
				mc_blas_vector_at(d, i) = mc_blas_vector_at(d, j);
				mc_blas_vector_at(d, j) = tmp;
				goto F60;
			}
			if (j - start > endd - j - 1) {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
			} else {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
			}
		} else {
			i = start - 1;
			j = endd + 1;
F90:
F100:
			j = j - 1;
			if (mc_blas_vector_at(d, j) > dmnmx) {
				goto F100;
			}
	F110:
			i = i + 1;
			if (mc_blas_vector_at(d, i) < dmnmx) {
				goto F110;
			}
			if (i < j) {
				tmp                     = mc_blas_vector_at(d, i);
				mc_blas_vector_at(d, i) = mc_blas_vector_at(d, j);
				mc_blas_vector_at(d, j) = tmp;
				goto F90;
			}
			if (j - start > endd - j - 1) {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
			} else {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
			}
		}
	}
	if (stkpnt > 0) {
		goto F10;
	}
}

#pragma mark - mc_lapack_llasrt -

MC_TARGET_FUNC void mc_lapack_llasrt(const char id, const int n, long double * d, int * info)
{
	int select = 20;
	int dir, endd, i, j, start, stkpnt;
	long double d1, d2, d3, dmnmx, tmp;

	int stack[2 * 32] = { 0 };

	*info = 0;
	 dir  = -1;
	if (mc_blas_lsame(id, 'D')) {
		dir = 0;
	} else if (mc_blas_lsame(id, 'I')) {
		dir = 1;
	}
	if (dir == -1) {
		*info = -1;
	} else if (n < 0) {
		*info = -2;
	}
	if (*info != 0) {
		mc_blas_xerbla("LLASRT", -(*info));
		return;
	}

	if (n <= 1) {
		return;
	}

	stkpnt                                = 1;
	mc_blas_matrix_at(stack, 2, 32, 1, 1) = 1;
	mc_blas_matrix_at(stack, 2, 32, 2, 1) = n;

F10:
	start  = mc_blas_matrix_at(stack, 2, 32, 1, stkpnt);
	endd   = mc_blas_matrix_at(stack, 2, 32, 2, stkpnt);
	stkpnt = stkpnt - 1;
	if (endd - start <= select && endd - start > 0) {
		if (dir == 0) {
			for (i = (start + 1); i <= endd; ++i) {
				for (j = i; j >= (start + 1); --j) {
					if (mc_blas_vector_at(d, j) > mc_blas_vector_at(d, j - 1)) {
						dmnmx                       = mc_blas_vector_at(d, j);
						mc_blas_vector_at(d, j)     = mc_blas_vector_at(d, j - 1);
						mc_blas_vector_at(d, j - 1) = dmnmx;
					} else {
						goto F30;
					}
				}
F30:;
			}
		} else {
			for (i = (start + 1); i <= endd; ++i) {
				for (j = i; j >= (start + 1); --j) {
					if (mc_blas_vector_at(d, j) < mc_blas_vector_at(d, j - 1)) {
						dmnmx                       = mc_blas_vector_at(d, j);
						mc_blas_vector_at(d, j)     = mc_blas_vector_at(d, j - 1);
						mc_blas_vector_at(d, j - 1) = dmnmx;
					} else {
						goto F50;
					}
				}
F50:;
			}
		}

	} else if (endd - start > select) {
		d1 = mc_blas_vector_at(d, start);
		d2 = mc_blas_vector_at(d, endd);
		i  = (start + endd) / 2;
		d3 = mc_blas_vector_at(d, i);
		if (d1 < d2) {
			if (d3 < d1) {
				dmnmx = d1;
			} else if (d3 < d2) {
				dmnmx = d3;
			} else {
				dmnmx = d2;
			}
		} else {
			if (d3 < d2) {
				dmnmx = d2;
			} else if (d3 < d1) {
				dmnmx = d3;
			} else {
				dmnmx = d1;
			}
		}
		if (dir == 0) {
			i = start - 1;
			j = endd + 1;
F60:
F70:
			j = j - 1;
			if (mc_blas_vector_at(d, j) < dmnmx) {
				goto F70;
			}
F80:
			i = i + 1;
			if (mc_blas_vector_at(d, i) > dmnmx) {
				goto F80;
			}
			if (i < j) {
				tmp                     = mc_blas_vector_at(d, i);
				mc_blas_vector_at(d, i) = mc_blas_vector_at(d, j);
				mc_blas_vector_at(d, j) = tmp;
				goto F60;
			}
			if (j - start > endd - j - 1) {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
			} else {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
			}
		} else {
			i = start - 1;
			j = endd + 1;
F90:
F100:
			j = j - 1;
			if (mc_blas_vector_at(d, j) > dmnmx) {
				goto F100;
			}
	F110:
			i = i + 1;
			if (mc_blas_vector_at(d, i) < dmnmx) {
				goto F110;
			}
			if (i < j) {
				tmp                     = mc_blas_vector_at(d, i);
				mc_blas_vector_at(d, i) = mc_blas_vector_at(d, j);
				mc_blas_vector_at(d, j) = tmp;
				goto F90;
			}
			if (j - start > endd - j - 1) {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
			} else {
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = j + 1;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = endd;
				stkpnt                                     = stkpnt + 1;
				mc_blas_matrix_at(stack, 2, 32, 1, stkpnt) = start;
				mc_blas_matrix_at(stack, 2, 32, 2, stkpnt) = j;
			}
		}
	}
	if (stkpnt > 0) {
		goto F10;
	}
}

#endif /* !MC_LAPACKE_LASRT_H */

/* EOF */