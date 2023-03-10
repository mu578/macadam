//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasdt.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>
#include <macadam/details/math/mc_itrunc32.h>
#include <macadam/details/math/mc_log.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_LAPACKE_LASDT_H
#define MC_LAPACKE_LASDT_H

#pragma mark - mc_lapack_slasdt -

MC_TARGET_FUNC void mc_lapack_slasdt(const int n, int * lvl, int * nd, int * inode, int * ndiml, int * ndimr, const int msub)
{
	const float two = 2.0f;

	int i, il, ir, llst, maxn, ncrnt, nlvl;
	float temp;

	 maxn                       = mc_maxmag(1, n);
	 temp                       = mc_logf(mc_cast(float, maxn) / mc_cast(const float, (msub + 1))) / mc_logf(two);
	*lvl                        = mc_itrunc32f(temp) + 1;

	 i                          = n / 2;
	mc_blas_vector_at(inode, 1) = i + 1;
	mc_blas_vector_at(ndiml, 1) = i;
	mc_blas_vector_at(ndimr, 1) = n - i - 1;
	 il                         = 0;
	 ir                         = 1;
	 llst                       = 1;
	for (nlvl = 1; nlvl <= (*lvl - 1); ++nlvl) {
		for (i = 0; i <= (llst - 1); ++i) {
			il                           = il + 2;
			ir                           = ir + 2;
			ncrnt                        = llst + i;
			mc_blas_vector_at(ndiml, il) = mc_blas_vector_at(ndiml, ncrnt) / 2;
			mc_blas_vector_at(ndimr, il) = mc_blas_vector_at(ndiml, ncrnt) - mc_blas_vector_at(ndiml, il) - 1;
			mc_blas_vector_at(inode, il) = mc_blas_vector_at(inode, ncrnt) - mc_blas_vector_at(ndimr, il) - 1;
			mc_blas_vector_at(ndiml, ir) = mc_blas_vector_at(ndimr, ncrnt) / 2;
			mc_blas_vector_at(ndimr, ir) = mc_blas_vector_at(ndimr, ncrnt) - mc_blas_vector_at(ndiml, ir) - 1;
			mc_blas_vector_at(inode, ir) = mc_blas_vector_at(inode, ncrnt) + mc_blas_vector_at(ndiml, ir) + 1;
		}
		llst = llst * 2;
	}
	*nd = (llst * 2) - 1;
}

#pragma mark - mc_lapack_dlasdt -

MC_TARGET_FUNC void mc_lapack_dlasdt(const int n, int * lvl, int * nd, int * inode, int * ndiml, int * ndimr, const int msub)
{
	const double two = 2.0;

	int i, il, ir, llst, maxn, ncrnt, nlvl;
	double temp;

	 maxn                       = mc_maxmag(1, n);
	 temp                       = mc_log(mc_cast(double, maxn) / mc_cast(const double, (msub + 1))) / mc_log(two);
	*lvl                        = mc_itrunc32(temp) + 1;

	 i                          = n / 2;
	mc_blas_vector_at(inode, 1) = i + 1;
	mc_blas_vector_at(ndiml, 1) = i;
	mc_blas_vector_at(ndimr, 1) = n - i - 1;
	 il                         = 0;
	 ir                         = 1;
	 llst                       = 1;
	for (nlvl = 1; nlvl <= (*lvl - 1); ++nlvl) {
		for (i = 0; i <= (llst - 1); ++i) {
			il                           = il + 2;
			ir                           = ir + 2;
			ncrnt                        = llst + i;
			mc_blas_vector_at(ndiml, il) = mc_blas_vector_at(ndiml, ncrnt) / 2;
			mc_blas_vector_at(ndimr, il) = mc_blas_vector_at(ndiml, ncrnt) - mc_blas_vector_at(ndiml, il) - 1;
			mc_blas_vector_at(inode, il) = mc_blas_vector_at(inode, ncrnt) - mc_blas_vector_at(ndimr, il) - 1;
			mc_blas_vector_at(ndiml, ir) = mc_blas_vector_at(ndimr, ncrnt) / 2;
			mc_blas_vector_at(ndimr, ir) = mc_blas_vector_at(ndimr, ncrnt) - mc_blas_vector_at(ndiml, ir) - 1;
			mc_blas_vector_at(inode, ir) = mc_blas_vector_at(inode, ncrnt) + mc_blas_vector_at(ndiml, ir) + 1;
		}
		llst = llst * 2;
	}
	*nd = (llst * 2) - 1;
}

#pragma mark - mc_lapack_llasdt -

MC_TARGET_FUNC void mc_lapack_llasdt(const int n, int * lvl, int * nd, int * inode, int * ndiml, int * ndimr, const int msub)
{
	const long double two = 2.0L;

	int i, il, ir, llst, maxn, ncrnt, nlvl;
	long double temp;

	 maxn                       = mc_maxmag(1, n);
	 temp                       = mc_logl(mc_cast(long double, maxn) / mc_cast(const long double, (msub + 1))) / mc_logl(two);
	*lvl                        = mc_itrunc32l(temp) + 1;

	 i                          = n / 2;
	mc_blas_vector_at(inode, 1) = i + 1;
	mc_blas_vector_at(ndiml, 1) = i;
	mc_blas_vector_at(ndimr, 1) = n - i - 1;
	 il                         = 0;
	 ir                         = 1;
	 llst                       = 1;
	for (nlvl = 1; nlvl <= (*lvl - 1); ++nlvl) {
		for (i = 0; i <= (llst - 1); ++i) {
			il                           = il + 2;
			ir                           = ir + 2;
			ncrnt                        = llst + i;
			mc_blas_vector_at(ndiml, il) = mc_blas_vector_at(ndiml, ncrnt) / 2;
			mc_blas_vector_at(ndimr, il) = mc_blas_vector_at(ndiml, ncrnt) - mc_blas_vector_at(ndiml, il) - 1;
			mc_blas_vector_at(inode, il) = mc_blas_vector_at(inode, ncrnt) - mc_blas_vector_at(ndimr, il) - 1;
			mc_blas_vector_at(ndiml, ir) = mc_blas_vector_at(ndimr, ncrnt) / 2;
			mc_blas_vector_at(ndimr, ir) = mc_blas_vector_at(ndimr, ncrnt) - mc_blas_vector_at(ndiml, ir) - 1;
			mc_blas_vector_at(inode, ir) = mc_blas_vector_at(inode, ncrnt) + mc_blas_vector_at(ndiml, ir) + 1;
		}
		llst = llst * 2;
	}
	*nd = (llst * 2) - 1;
}

#endif /* !MC_LAPACKE_LASDT_H */

/* EOF */