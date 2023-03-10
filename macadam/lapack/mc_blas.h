//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#ifndef MC_BLAS_H
#define MC_BLAS_H

#	undef  MC_TARGET_BLAS_USE_CLAYOUT
#	undef  MC_TARGET_BLAS_USE_FLAYOUT
#	define MC_TARGET_BLAS_USE_CLAYOUT 1
#	define MC_TARGET_BLAS_USE_FLAYOUT 0

#	include <macadam/lapack/blas/mc_blas_abs1.h>
#	include <macadam/lapack/blas/mc_blas_access.h>
#	include <macadam/lapack/blas/mc_blas_asum.h>
#	include <macadam/lapack/blas/mc_blas_axpy.h>
#	include <macadam/lapack/blas/mc_blas_copy.h>
#	include <macadam/lapack/blas/mc_blas_dot.h>
#	include <macadam/lapack/blas/mc_blas_gbmv.h>
#	include <macadam/lapack/blas/mc_blas_gemm.h>
#	include <macadam/lapack/blas/mc_blas_gemv.h>
#	include <macadam/lapack/blas/mc_blas_ger.h>
#	include <macadam/lapack/blas/mc_blas_iamax.h>
#	include <macadam/lapack/blas/mc_blas_lsame.h>
#	include <macadam/lapack/blas/mc_blas_nrm2.h>
#	include <macadam/lapack/blas/mc_blas_rot.h>
#	include <macadam/lapack/blas/mc_blas_rotg.h>
#	include <macadam/lapack/blas/mc_blas_rotm.h>
#	include <macadam/lapack/blas/mc_blas_rotmg.h>
#	include <macadam/lapack/blas/mc_blas_sbmv.h>
#	include <macadam/lapack/blas/mc_blas_scal.h>
#	include <macadam/lapack/blas/mc_blas_spmv.h>
#	include <macadam/lapack/blas/mc_blas_spr.h>
#	include <macadam/lapack/blas/mc_blas_spr2.h>
#	include <macadam/lapack/blas/mc_blas_symm.h>
#	include <macadam/lapack/blas/mc_blas_symv.h>
#	include <macadam/lapack/blas/mc_blas_syr.h>
#	include <macadam/lapack/blas/mc_blas_syr2.h>
#	include <macadam/lapack/blas/mc_blas_syr2k.h>
#	include <macadam/lapack/blas/mc_blas_syrk.h>
#	include <macadam/lapack/blas/mc_blas_tbmv.h>
#	include <macadam/lapack/blas/mc_blas_tbsv.h>
#	include <macadam/lapack/blas/mc_blas_tpmv.h>
#	include <macadam/lapack/blas/mc_blas_tpsv.h>
#	include <macadam/lapack/blas/mc_blas_trmm.h>
#	include <macadam/lapack/blas/mc_blas_trmv.h>
#	include <macadam/lapack/blas/mc_blas_trsm.h>
#	include <macadam/lapack/blas/mc_blas_trsv.h>
#	include <macadam/lapack/blas/mc_blas_xerbla.h>

#	if MC_TARGET_BLAS_USE_NATIVE

#		include <macadam/lapack/blas/native/mc_blas_native_abs1.h>
#		include <macadam/lapack/blas/native/mc_blas_native_asum.h>
#		include <macadam/lapack/blas/native/mc_blas_native_axpy.h>
#		include <macadam/lapack/blas/native/mc_blas_native_copy.h>
#		include <macadam/lapack/blas/native/mc_blas_native_dot.h>
#		include <macadam/lapack/blas/native/mc_blas_native_gbmv.h>
#		include <macadam/lapack/blas/native/mc_blas_native_gemm.h>
#		include <macadam/lapack/blas/native/mc_blas_native_gemv.h>
#		include <macadam/lapack/blas/native/mc_blas_native_ger.h>
#		include <macadam/lapack/blas/native/mc_blas_native_iamax.h>
#		include <macadam/lapack/blas/native/mc_blas_native_nrm2.h>
#		include <macadam/lapack/blas/native/mc_blas_native_rot.h>
#		include <macadam/lapack/blas/native/mc_blas_native_rotg.h>
#		include <macadam/lapack/blas/native/mc_blas_native_rotm.h>
#		include <macadam/lapack/blas/native/mc_blas_native_rotmg.h>
#		include <macadam/lapack/blas/native/mc_blas_native_sbmv.h>
#		include <macadam/lapack/blas/native/mc_blas_native_scal.h>
#		include <macadam/lapack/blas/native/mc_blas_native_spmv.h>
#		include <macadam/lapack/blas/native/mc_blas_native_spr.h>
#		include <macadam/lapack/blas/native/mc_blas_native_spr2.h>
#		include <macadam/lapack/blas/native/mc_blas_native_swap.h>
#		include <macadam/lapack/blas/native/mc_blas_native_symm.h>
#		include <macadam/lapack/blas/native/mc_blas_native_symv.h>
#		include <macadam/lapack/blas/native/mc_blas_native_syr.h>
#		include <macadam/lapack/blas/native/mc_blas_native_syr2.h>
#		include <macadam/lapack/blas/native/mc_blas_native_syr2k.h>
#		include <macadam/lapack/blas/native/mc_blas_native_syrk.h>
#		include <macadam/lapack/blas/native/mc_blas_native_tbmv.h>

#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		define mc_scabs1 mc_blas_native_scabs1
#		define mc_dcabs1 mc_blas_native_dcabs1
#		define mc_szabs1 mc_blas_native_szabs1
#		define mc_dzabs1 mc_blas_native_dzabs1
#	else
#		define mc_scabs1 mc_blas_scabs1
#		define mc_dcabs1 mc_blas_dcabs1
#		define mc_szabs1 mc_blas_szabs1
#		define mc_dzabs1 mc_blas_dzabs1
#	endif
#		define mc_sasum  mc_blas_native_sasum
#		define mc_dasum  mc_blas_native_dasum
#		define mc_scasum mc_blas_native_scasum
#		define mc_dzasum mc_blas_native_dzasum
#		define mc_saxpy  mc_blas_native_saxpy
#		define mc_daxpy  mc_blas_native_daxpy
#		define mc_caxpy  mc_blas_native_caxpy
#		define mc_zaxpy  mc_blas_native_zaxpy
#		define mc_scopy  mc_blas_native_scopy
#		define mc_dcopy  mc_blas_native_dcopy
#		define mc_ccopy  mc_blas_native_ccopy
#		define mc_zcopy  mc_blas_native_zcopy
#		define mc_sdot   mc_blas_native_sdot
#		define mc_ddot   mc_blas_native_ddot
#		define mc_cdotc  mc_blas_native_cdotc
#		define mc_zdotc  mc_blas_native_zdotc
#		define mc_cdotu  mc_blas_native_cdotu
#		define mc_zdotu  mc_blas_native_zdotu
#		define mc_sgbmv  mc_blas_native_sgbmv
#		define mc_dgbmv  mc_blas_native_dgbmv
#		define mc_cgbmv  mc_blas_native_cgbmv
#		define mc_zgbmv  mc_blas_native_zgbmv
#		define mc_sgemm  mc_blas_native_sgemm
#		define mc_dgemm  mc_blas_native_dgemm
#		define mc_cgemm  mc_blas_native_cgemm
#		define mc_zgemm  mc_blas_native_zgemm
#		define mc_sgemv  mc_blas_native_sgemv
#		define mc_dgemv  mc_blas_native_dgemv
#		define mc_cgemv  mc_blas_native_cgemv
#		define mc_zgemv  mc_blas_native_zgemv
#		define mc_sger   mc_blas_native_sger
#		define mc_dger   mc_blas_native_dger
#		define mc_cgerc  mc_blas_native_cgerc
#		define mc_zgerc  mc_blas_native_zgerc
#		define mc_cgeru  mc_blas_native_cgeru
#		define mc_zgeru  mc_blas_native_zgeru
#		define mc_isamax mc_blas_native_isamax
#		define mc_idamax mc_blas_native_idamax
#		define mc_icamax mc_blas_native_icamax
#		define mc_izamax mc_blas_native_izamax
#		define mc_snrm2  mc_blas_native_snrm2
#		define mc_dnrm2  mc_blas_native_dnrm2
#		define mc_scnrm2 mc_blas_native_scnrm2
#		define mc_dznrm2 mc_blas_native_dznrm2
#		define mc_srot   mc_blas_native_srot
#		define mc_drot   mc_blas_native_drot
#	if !MC_TARGET_BLAS_USE_OPENBLAS
#		define mc_csrot  mc_blas_native_csrot
#		define mc_zdrot  mc_blas_native_zdrot
#	else
#		define mc_csrot  mc_blas_csrot
#		define mc_zdrot  mc_blas_zdrot
#	endif
#		define mc_srotg  mc_blas_native_srotg
#		define mc_drotg  mc_blas_native_drotg
#	if !MC_TARGET_BLAS_USE_OPENBLAS
#		define mc_crotg  mc_blas_native_crotg
#		define mc_zrotg  mc_blas_native_zrotg
#	else
#		define mc_crotg  mc_blas_crotg
#		define mc_zrotg  mc_blas_zrotg
#	endif
#		define mc_srotm  mc_blas_native_srotm
#		define mc_drotm  mc_blas_native_drotm
#		define mc_srotmg mc_blas_native_srotmg
#		define mc_drotmg mc_blas_native_drotmg
#		define mc_ssbmv  mc_blas_native_ssbmv
#		define mc_dsbmv  mc_blas_native_dsbmv
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		define mc_csbmv  mc_blas_native_csbmv
#		define mc_zsbmv  mc_blas_native_zsbmv
#	else
#		define mc_csbmv  mc_blas_csbmv
#		define mc_zsbmv  mc_blas_zsbmv
#	endif
#		define mc_sscal  mc_blas_native_sscal
#		define mc_dscal  mc_blas_native_dscal
#		define mc_cscal  mc_blas_native_cscal
#		define mc_zscal  mc_blas_native_zscal
#		define mc_sspmv  mc_blas_native_sspmv
#		define mc_dspmv  mc_blas_native_dspmv
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		define mc_cspmv  mc_blas_native_cspmv
#		define mc_zspmv  mc_blas_native_zspmv
#	else
#		define mc_cspmv  mc_blas_cspmv
#		define mc_zspmv  mc_blas_zspmv
#	endif
#		define mc_sspr   mc_blas_native_sspr
#		define mc_dspr   mc_blas_native_dspr
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		define mc_cspr   mc_blas_native_cspr
#		define mc_zspr   mc_blas_native_zspr
#	else
#		define mc_cspr   mc_blas_cspr
#		define mc_zspr   mc_blas_zspr
#	endif
#		define mc_sspr2  mc_blas_native_sspr2
#		define mc_dspr2  mc_blas_native_dspr2
#		define mc_sswap  mc_blas_native_sswap
#		define mc_dswap  mc_blas_native_dswap
#		define mc_cswap  mc_blas_native_cswap
#		define mc_zswap  mc_blas_native_zswap
#		define mc_ssymm  mc_blas_native_ssymm
#		define mc_dsymm  mc_blas_native_dsymm
#		define mc_csymm  mc_blas_native_csymm
#		define mc_zsymm  mc_blas_native_zsymm
#		define mc_ssymv  mc_blas_native_ssymv
#		define mc_dsymv  mc_blas_native_dsymv
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		define mc_csymv  mc_blas_native_csymv
#		define mc_zsymv  mc_blas_native_zsymv
#	else
#		define mc_csymv  mc_blas_csymv
#		define mc_zsymv  mc_blas_zsymv
#	endif
#		define mc_ssyr   mc_blas_native_ssyr
#		define mc_dsyr   mc_blas_native_dsyr
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		define mc_csyr   mc_blas_native_csyr
#		define mc_zsyr   mc_blas_native_zsyr
#	else
#		define mc_csyr   mc_blas_csyr
#		define mc_zsyr   mc_blas_zsyr
#	endif
#		define mc_ssyr2  mc_blas_native_ssyr2
#		define mc_dsyr2  mc_blas_native_dsyr2
#	if !MC_TARGET_BLAS_USE_OPENBLAS   \
	&& !MC_TARGET_BLAS_USE_ACCELERATE \
	&& !MC_TARGET_BLAS_USE_VECLIB
#		define mc_csyr2  mc_blas_native_csyr2
#		define mc_zsyr2  mc_blas_native_zsyr2
#	else
#		define mc_csyr2  mc_blas_csyr2
#		define mc_zsyr2  mc_blas_zsyr2
#	endif
#		define mc_ssyr2k mc_blas_native_ssyr2k
#		define mc_dsyr2k mc_blas_native_dsyr2k
#		define mc_csyr2k mc_blas_native_csyr2k
#		define mc_zsyr2k mc_blas_native_zsyr2k
#		define mc_ssyrk  mc_blas_native_ssyrk
#		define mc_dsyrk  mc_blas_native_dsyrk
#		define mc_csyrk  mc_blas_native_csyrk
#		define mc_zsyrk  mc_blas_native_zsyrk
#		define mc_stbmv  mc_blas_native_stbmv
#		define mc_dtbmv  mc_blas_native_dtbmv
#		define mc_ctbmv  mc_blas_native_ctbmv
#		define mc_ztbmv  mc_blas_native_ztbmv

#	else

#		define mc_scabs1 mc_blas_scabs1
#		define mc_dcabs1 mc_blas_dcabs1
#		define mc_szabs1 mc_blas_szabs1
#		define mc_dzabs1 mc_blas_dzabs1
#		define mc_sasum  mc_blas_sasum
#		define mc_dasum  mc_blas_dasum
#		define mc_scasum mc_blas_scasum
#		define mc_dzasum mc_blas_dzasum
#		define mc_saxpy  mc_blas_saxpy
#		define mc_daxpy  mc_blas_daxpy
#		define mc_caxpy  mc_blas_caxpy
#		define mc_zaxpy  mc_blas_zaxpy
#		define mc_scopy  mc_blas_scopy
#		define mc_dcopy  mc_blas_dcopy
#		define mc_ccopy  mc_blas_ccopy
#		define mc_zcopy  mc_blas_zcopy
#		define mc_sdot   mc_blas_sdot
#		define mc_ddot   mc_blas_ddot
#		define mc_cdotc  mc_blas_cdotc
#		define mc_zdotc  mc_blas_zdotc
#		define mc_cdotu  mc_blas_cdotu
#		define mc_zdotu  mc_blas_zdotu
#		define mc_sgbmv  mc_blas_sgbmv
#		define mc_dgbmv  mc_blas_dgbmv
#		define mc_cgbmv  mc_blas_cgbmv
#		define mc_zgbmv  mc_blas_zgbmv
#		define mc_sgemm  mc_blas_sgemm
#		define mc_dgemm  mc_blas_dgemm
#		define mc_cgemm  mc_blas_cgemm
#		define mc_zgemm  mc_blas_zgemm
#		define mc_sgemv  mc_blas_sgemv
#		define mc_dgemv  mc_blas_dgemv
#		define mc_cgemv  mc_blas_cgemv
#		define mc_zgemv  mc_blas_zgemv
#		define mc_sger   mc_blas_sger
#		define mc_dger   mc_blas_dger
#		define mc_cgerc  mc_blas_cgerc
#		define mc_zgerc  mc_blas_zgerc
#		define mc_cgeru  mc_blas_cgeru
#		define mc_zgeru  mc_blas_zgeru
#		define mc_isamax mc_blas_isamax
#		define mc_idamax mc_blas_idamax
#		define mc_icamax mc_blas_icamax
#		define mc_izamax mc_blas_izamax
#		define mc_snrm2  mc_blas_snrm2
#		define mc_dnrm2  mc_blas_dnrm2
#		define mc_scnrm2 mc_blas_scnrm2
#		define mc_dznrm2 mc_blas_dznrm2
#		define mc_srot   mc_blas_srot
#		define mc_drot   mc_blas_drot
#		define mc_csrot  mc_blas_csrot
#		define mc_zdrot  mc_blas_zdrot
#		define mc_srotg  mc_blas_srotg
#		define mc_drotg  mc_blas_drotg
#		define mc_crotg  mc_blas_crotg
#		define mc_zrotg  mc_blas_zrotg
#		define mc_srotm  mc_blas_srotm
#		define mc_drotm  mc_blas_drotm
#		define mc_srotmg mc_blas_srotmg
#		define mc_drotmg mc_blas_drotmg
#		define mc_ssbmv  mc_blas_ssbmv
#		define mc_dsbmv  mc_blas_dsbmv
#		define mc_csbmv  mc_blas_csbmv
#		define mc_zsbmv  mc_blas_zsbmv
#		define mc_sscal  mc_blas_sscal
#		define mc_dscal  mc_blas_dscal
#		define mc_cscal  mc_blas_cscal
#		define mc_zscal  mc_blas_zscal
#		define mc_sspmv  mc_blas_sspmv
#		define mc_dspmv  mc_blas_dspmv
#		define mc_cspmv  mc_blas_cspmv
#		define mc_zspmv  mc_blas_zspmv
#		define mc_sspr   mc_blas_sspr
#		define mc_dspr   mc_blas_dspr
#		define mc_cspr   mc_blas_cspr
#		define mc_zspr   mc_blas_zspr
#		define mc_sspr2  mc_blas_sspr2
#		define mc_dspr2  mc_blas_dspr2
#		define mc_sswap  mc_blas_sswap
#		define mc_dswap  mc_blas_dswap
#		define mc_cswap  mc_blas_cswap
#		define mc_zswap  mc_blas_zswap
#		define mc_ssymm  mc_blas_ssymm
#		define mc_dsymm  mc_blas_dsymm
#		define mc_csymm  mc_blas_csymm
#		define mc_zsymm  mc_blas_zsymm
#		define mc_ssymv  mc_blas_ssymv
#		define mc_dsymv  mc_blas_dsymv
#		define mc_csymv  mc_blas_csymv
#		define mc_zsymv  mc_blas_zsymv
#		define mc_ssyr   mc_blas_ssyr
#		define mc_dsyr   mc_blas_dsyr
#		define mc_csyr   mc_blas_csyr
#		define mc_zsyr   mc_blas_zsyr
#		define mc_ssyr2  mc_blas_ssyr2
#		define mc_dsyr2  mc_blas_dsyr2
#		define mc_csyr2  mc_blas_csyr2
#		define mc_zsyr2  mc_blas_zsyr2
#		define mc_ssyr2k mc_blas_ssyr2k
#		define mc_dsyr2k mc_blas_dsyr2k
#		define mc_csyr2k mc_blas_csyr2k
#		define mc_zsyr2k mc_blas_zsyr2k
#		define mc_ssyrk  mc_blas_ssyrk
#		define mc_dsyrk  mc_blas_dsyrk
#		define mc_csyrk  mc_blas_csyrk
#		define mc_zsyrk  mc_blas_zsyrk
#		define mc_stbmv  mc_blas_stbmv
#		define mc_dtbmv  mc_blas_dtbmv
#		define mc_ctbmv  mc_blas_ctbmv
#		define mc_ztbmv  mc_blas_ztbmv

#	endif

#	define mc_sdsasum mc_blas_sdsasum
#	define mc_lasum   mc_blas_lasum
#	define mc_lqasum  mc_blas_lqasum
#	define mc_laxpy   mc_blas_laxpy
#	define mc_qaxpy   mc_blas_qaxpy
#	define mc_lcopy   mc_blas_lcopy
#	define mc_qcopy   mc_blas_qcopy
#	define mc_ldot    mc_blas_ldot
#	define mc_qdotc   mc_blas_qdotc
#	define mc_qdotu   mc_blas_qdotu
#	define mc_lgbmv   mc_blas_lgbmv
#	define mc_qgbmv   mc_blas_qgbmv
#	define mc_lgemm   mc_blas_lgemm
#	define mc_qgemm   mc_blas_qgemm
#	define mc_lgemv   mc_blas_lgemv
#	define mc_qgemv   mc_blas_qgemv
#	define mc_lger    mc_blas_lger
#	define mc_qgerc   mc_blas_qgerc
#	define mc_qgeru   mc_blas_qgeru
#	define mc_ilamax  mc_blas_ilamax
#	define mc_iqamax  mc_blas_iqamax
#	define mc_dsnrm2  mc_blas_dsnrm2
#	define mc_lnrm2   mc_blas_lnrm2
#	define mc_lqnrm2  mc_blas_lqnrm2
#	define mc_lrot    mc_blas_lrot
#	define mc_qlrot   mc_blas_qlrot
#	define mc_lrotg   mc_blas_lrotg
#	define mc_qrotg   mc_blas_qrotg
#	define mc_lrotm   mc_blas_lrotm
#	define mc_lrotmg  mc_blas_lrotmg
#	define mc_lsbmv   mc_blas_lsbmv
#	define mc_qsbmv   mc_blas_qsbmv
#	define mc_lscal   mc_blas_lscal
#	define mc_qscal   mc_blas_qscal
#	define mc_lspmv   mc_blas_lspmv
#	define mc_qspmv   mc_blas_qspmv
#	define mc_lspr    mc_blas_lspr
#	define mc_qspr    mc_blas_qspr
#	define mc_lspr2   mc_blas_lspr2
#	define mc_lswap   mc_blas_lswap
#	define mc_qswap   mc_blas_qswap
#	define mc_lsymm   mc_blas_lsymm
#	define mc_qsymm   mc_blas_qsymm
#	define mc_lsymv   mc_blas_lsymv
#	define mc_qsymv   mc_blas_qsymv
#	define mc_lsyr    mc_blas_lsyr
#	define mc_qsyr    mc_blas_qsyr
#	define mc_lsyr2   mc_blas_lsyr2
#	define mc_qsyr2   mc_blas_qsyr2
#	define mc_lsyr2k  mc_blas_lsyr2k
#	define mc_qsyr2k  mc_blas_qsyr2k
#	define mc_lsyrk   mc_blas_lsyrk
#	define mc_qsyrk   mc_blas_qsyrk
#	define mc_ltbmv   mc_blas_ltbmv
#	define mc_qtbmv   mc_blas_qtbmv

#endif /* !MC_BLAS_H */

/* EOF */