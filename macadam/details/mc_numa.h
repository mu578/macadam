//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_numa.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#ifndef MC_NUMA_H
#define MC_NUMA_H

#	include <macadam/details/numa/mc_2summx1.h>
#	include <macadam/details/numa/mc_2sum1xn.h>
#	include <macadam/details/numa/mc_a2summx1.h>
#	include <macadam/details/numa/mc_a2sum1xn.h>
#	include <macadam/details/numa/mc_addmulatx2x2.h>
#	include <macadam/details/numa/mc_addmulatx3x3.h>
#	include <macadam/details/numa/mc_addmulax2x2.h>
#	include <macadam/details/numa/mc_addmulax3x3.h>
#	include <macadam/details/numa/mc_addxv1xn.h>
#	include <macadam/details/numa/mc_addxy1xn.h>
#	include <macadam/details/numa/mc_amean1xn.h>
#	include <macadam/details/numa/mc_arange1xn.h>
#	include <macadam/details/numa/mc_argmax1xn.h>
#	include <macadam/details/numa/mc_argmaxmx1.h>
#	include <macadam/details/numa/mc_argmin1xn.h>
#	include <macadam/details/numa/mc_argminmx1.h>
#	include <macadam/details/numa/mc_argsort1xn.h>
#	include <macadam/details/numa/mc_asum1xn.h>
#	include <macadam/details/numa/mc_centermxn.h>
#	include <macadam/details/numa/mc_centernxn.h>
#	include <macadam/details/numa/mc_chpoly2x2.h>
#	include <macadam/details/numa/mc_chpoly3x3.h>
#	include <macadam/details/numa/mc_copy1xn.h>
#	include <macadam/details/numa/mc_copy2x2.h>
#	include <macadam/details/numa/mc_copy3x3.h>
#	include <macadam/details/numa/mc_copy4x4.h>
#	include <macadam/details/numa/mc_copymxn.h>
#	include <macadam/details/numa/mc_copynxn.h>
#	include <macadam/details/numa/mc_covar1xn.h>
#	include <macadam/details/numa/mc_covarmxn.h>
#	include <macadam/details/numa/mc_cumprod1xn.h>
#	include <macadam/details/numa/mc_cumsum1xn.h>
#	include <macadam/details/numa/mc_det2x2.h>
#	include <macadam/details/numa/mc_det3x3.h>
#	include <macadam/details/numa/mc_det4x4.h>
#	include <macadam/details/numa/mc_diag1xn.h>
#	include <macadam/details/numa/mc_diff1xn.h>
#	include <macadam/details/numa/mc_diffk1xn.h>
#	include <macadam/details/numa/mc_divxv1xn.h>
#	include <macadam/details/numa/mc_divxy1xn.h>
#	include <macadam/details/numa/mc_dotp1xn.h>
#	include <macadam/details/numa/mc_dotp3x1.h>
#	include <macadam/details/numa/mc_dotpmx1.h>
#	include <macadam/details/numa/mc_dotpmxn.h>
#	include <macadam/details/numa/mc_dotpnxn.h>
#	include <macadam/details/numa/mc_eig2x2.h>
#	include <macadam/details/numa/mc_eigsy2x2.h>
#	include <macadam/details/numa/mc_eigsy3x3.h>
#	include <macadam/details/numa/mc_eye2x2.h>
#	include <macadam/details/numa/mc_eye3x3.h>
#	include <macadam/details/numa/mc_eyenxn.h>
#	include <macadam/details/numa/mc_flip1xn.h>
#	include <macadam/details/numa/mc_fliplrmxn.h>
#	include <macadam/details/numa/mc_flipmx1.h>
#	include <macadam/details/numa/mc_fliprg1xn.h>
#	include <macadam/details/numa/mc_flipudmxn.h>
#	include <macadam/details/numa/mc_frobnormmxn.h>
#	include <macadam/details/numa/mc_frobnormnxn.h>
#	include <macadam/details/numa/mc_histcs1xn.h>
#	include <macadam/details/numa/mc_histce1xn.h>
#	include <macadam/details/numa/mc_histcg1xn.h>
#	include <macadam/details/numa/mc_histcb1xn.h>
#	include <macadam/details/numa/mc_icopy1xn.h>
#	include <macadam/details/numa/mc_infnorm1xn.h>
#	include <macadam/details/numa/mc_infnormmxn.h>
#	include <macadam/details/numa/mc_infnormnxn.h>
#	include <macadam/details/numa/mc_isort1xn.h>
#	include <macadam/details/numa/mc_izeros1xn.h>
#	include <macadam/details/numa/mc_jacobisynxn.h>
#	include <macadam/details/numa/mc_kronmxn.h>
#	include <macadam/details/numa/mc_l1norm1xn.h>
#	include <macadam/details/numa/mc_l2norm1x2.h>
#	include <macadam/details/numa/mc_l2norm1x3.h>
#	include <macadam/details/numa/mc_l2norm1xn.h>
#	include <macadam/details/numa/mc_l2norm2x1.h>
#	include <macadam/details/numa/mc_l2norm3x1.h>
#	include <macadam/details/numa/mc_l2normmx1.h>
#	include <macadam/details/numa/mc_ldu3x3.h>
#	include <macadam/details/numa/mc_ldup3x3.h>
#	include <macadam/details/numa/mc_linspace1xn.h>
#	include <macadam/details/numa/mc_logspace1xn.h>
#	include <macadam/details/numa/mc_lu2x2.h>
#	include <macadam/details/numa/mc_lu3x3.h>
#	include <macadam/details/numa/mc_lup3x3.h>
#	include <macadam/details/numa/mc_lupnxn.h>
#	include <macadam/details/numa/mc_lusolve3x3.h>
#	include <macadam/details/numa/mc_lusolvenxn.h>
#	include <macadam/details/numa/mc_magic2x2.h>
#	include <macadam/details/numa/mc_magic3x3.h>
#	include <macadam/details/numa/mc_magic4x4.h>
#	include <macadam/details/numa/mc_magic5x5.h>
#	include <macadam/details/numa/mc_magic6x6.h>
#	include <macadam/details/numa/mc_magicnxn.h>
#	include <macadam/details/numa/mc_max1xn.h>
#	include <macadam/details/numa/mc_maxmx1.h>
#	include <macadam/details/numa/mc_mean1xn.h>
#	include <macadam/details/numa/mc_meanmx1.h>
#	include <macadam/details/numa/mc_mgs3x3.h>
#	include <macadam/details/numa/mc_mgsmxn.h>
#	include <macadam/details/numa/mc_mgsnxn.h>
#	include <macadam/details/numa/mc_min1xn.h>
#	include <macadam/details/numa/mc_minmax1xn.h>
#	include <macadam/details/numa/mc_minmaxmx1.h>
#	include <macadam/details/numa/mc_minmx1.h>
#	include <macadam/details/numa/mc_minormxn.h>
#	include <macadam/details/numa/mc_moment1xn.h>
#	include <macadam/details/numa/mc_momentmx1.h>
#	include <macadam/details/numa/mc_mrmsmx1.h>
#	include <macadam/details/numa/mc_mssqr1xn.h>
#	include <macadam/details/numa/mc_mssqrmx1.h>
#	include <macadam/details/numa/mc_mstdd1xn.h>
#	include <macadam/details/numa/mc_mstddmx1.h>
#	include <macadam/details/numa/mc_mstddv1xn.h>
#	include <macadam/details/numa/mc_mulab2x2.h>
#	include <macadam/details/numa/mc_mulab3x3.h>
#	include <macadam/details/numa/mc_mulabmxn.h>
#	include <macadam/details/numa/mc_mulabnxn.h>
#	include <macadam/details/numa/mc_mulabt2x2.h>
#	include <macadam/details/numa/mc_mulabt3x3.h>
#	include <macadam/details/numa/mc_mulabtmxn.h>
#	include <macadam/details/numa/mc_mulabtnxn.h>
#	include <macadam/details/numa/mc_mulab2x2.h>
#	include <macadam/details/numa/mc_mulatb2x2.h>
#	include <macadam/details/numa/mc_mulatb3x3.h>
#	include <macadam/details/numa/mc_mulatbmxn.h>
#	include <macadam/details/numa/mc_mulatx2x2.h>
#	include <macadam/details/numa/mc_mulatx3x3.h>
#	include <macadam/details/numa/mc_mulatxmxn.h>
#	include <macadam/details/numa/mc_mulatxnxn.h>
#	include <macadam/details/numa/mc_mulax2x2.h>
#	include <macadam/details/numa/mc_mulax3x3.h>
#	include <macadam/details/numa/mc_mulaxmxn.h>
#	include <macadam/details/numa/mc_mulaxnxn.h>
#	include <macadam/details/numa/mc_muleabmxn.h>
#	include <macadam/details/numa/mc_muleabnxn.h>
#	include <macadam/details/numa/mc_mulxv1xn.h>
#	include <macadam/details/numa/mc_mulxy1xn.h>
#	include <macadam/details/numa/mc_nnz1xn.h>
#	include <macadam/details/numa/mc_nnzmx1.h>
#	include <macadam/details/numa/mc_ones1x2.h>
#	include <macadam/details/numa/mc_ones1x3.h>
#	include <macadam/details/numa/mc_ones1xn.h>
#	include <macadam/details/numa/mc_ones2x2.h>
#	include <macadam/details/numa/mc_ones3x3.h>
#	include <macadam/details/numa/mc_onesmx1.h>
#	include <macadam/details/numa/mc_onesmxn.h>
#	include <macadam/details/numa/mc_onesnxn.h>
#	include <macadam/details/numa/mc_ortho3x3.h>
#	include <macadam/details/numa/mc_orthomxn.h>
#	include <macadam/details/numa/mc_orthonxn.h>
#	include <macadam/details/numa/mc_orthrmxn.h>
#	include <macadam/details/numa/mc_orthrnxn.h>
#	include <macadam/details/numa/mc_outpabmxn.h>
#	include <macadam/details/numa/mc_outpxy2x2.h>
#	include <macadam/details/numa/mc_outpxy3x3.h>
#	include <macadam/details/numa/mc_outpxymxn.h>
#	include <macadam/details/numa/mc_outpxynxn.h>
#	include <macadam/details/numa/mc_poly2fit1xn.h>
#	include <macadam/details/numa/mc_poly3fit1xn.h>
#	include <macadam/details/numa/mc_poly4fit1xn.h>
#	include <macadam/details/numa/mc_poly5fit1xn.h>
#	include <macadam/details/numa/mc_poly6fit1xn.h>
#	include <macadam/details/numa/mc_poly7fit1xn.h>
#	include <macadam/details/numa/mc_poly8fit1xn.h>
#	include <macadam/details/numa/mc_poly9fit1xn.h>
#	include <macadam/details/numa/mc_polyfit1xn.h>
#	include <macadam/details/numa/mc_ppsdev1xn.h>
#	include <macadam/details/numa/mc_ppvar1xn.h>
#	include <macadam/details/numa/mc_qrgs2x2.h>
#	include <macadam/details/numa/mc_qrgs3x3.h>
#	include <macadam/details/numa/mc_qrgv3x3.h>
#	include <macadam/details/numa/mc_qrmgs3x3.h>
#	include <macadam/details/numa/mc_qrsolve3x3.h>
#	include <macadam/details/numa/mc_rescale1xn.h>
#	include <macadam/details/numa/mc_rescalemx1.h>
#	include <macadam/details/numa/mc_rescalemxn.h>
#	include <macadam/details/numa/mc_rescalenxn.h>
#	include <macadam/details/numa/mc_rmgsmxn.h>
#	include <macadam/details/numa/mc_rmgsnxn.h>
#	include <macadam/details/numa/mc_rms1xn.h>
#	include <macadam/details/numa/mc_rmsmx1.h>
#	include <macadam/details/numa/mc_rotl1xn.h>
#	include <macadam/details/numa/mc_sort1xn.h>
#	include <macadam/details/numa/mc_spsdev1xn.h>
#	include <macadam/details/numa/mc_spvar1xn.h>
#	include <macadam/details/numa/mc_ssqr1x2.h>
#	include <macadam/details/numa/mc_ssqr1x3.h>
#	include <macadam/details/numa/mc_ssqr1xn.h>
#	include <macadam/details/numa/mc_ssqr2x1.h>
#	include <macadam/details/numa/mc_ssqr3x1.h>
#	include <macadam/details/numa/mc_ssqrmx1.h>
#	include <macadam/details/numa/mc_stdd1xn.h>
#	include <macadam/details/numa/mc_stddmx1.h>
#	include <macadam/details/numa/mc_subxv1xn.h>
#	include <macadam/details/numa/mc_subxy1xn.h>
#	include <macadam/details/numa/mc_sum1xn.h>
#	include <macadam/details/numa/mc_summx1.h>
#	include <macadam/details/numa/mc_sumsq1xn.h>
#	include <macadam/details/numa/mc_svd2x2.h>
#	include <macadam/details/numa/mc_svd3x3.h>
#	include <macadam/details/numa/mc_svdmx3.h>
#	include <macadam/details/numa/mc_svdgr1mxn.h>
#	include <macadam/details/numa/mc_swap1xn.h>
#	include <macadam/details/numa/mc_swapmx1.h>
#	include <macadam/details/numa/mc_sytrize2x2.h>
#	include <macadam/details/numa/mc_sytrize3x3.h>
#	include <macadam/details/numa/mc_sytrizenxn.h>
#	include <macadam/details/numa/mc_trace2x2.h>
#	include <macadam/details/numa/mc_trace3x3.h>
#	include <macadam/details/numa/mc_tracemxn.h>
#	include <macadam/details/numa/mc_tracenxn.h>
#	include <macadam/details/numa/mc_trilnxn.h>
#	include <macadam/details/numa/mc_trilsolvenxn.h>
#	include <macadam/details/numa/mc_trilssqrnxn.h>
#	include <macadam/details/numa/mc_triunxn.h>
#	include <macadam/details/numa/mc_triusolve3x3.h>
#	include <macadam/details/numa/mc_triusolvenxn.h>
#	include <macadam/details/numa/mc_triussqrnxn.h>
#	include <macadam/details/numa/mc_trsi2x2.h>
#	include <macadam/details/numa/mc_trsi3x3.h>
#	include <macadam/details/numa/mc_trsimxn.h>
#	include <macadam/details/numa/mc_trsinxn.h>
#	include <macadam/details/numa/mc_trsp2x2.h>
#	include <macadam/details/numa/mc_trsp3x3.h>
#	include <macadam/details/numa/mc_trspmxn.h>
#	include <macadam/details/numa/mc_trspnxn.h>
#	include <macadam/details/numa/mc_ttvar1xn.h>
#	include <macadam/details/numa/mc_unit1x2.h>
#	include <macadam/details/numa/mc_unit1x3.h>
#	include <macadam/details/numa/mc_unit2x1.h>
#	include <macadam/details/numa/mc_unit3x1.h>
#	include <macadam/details/numa/mc_unitmxn.h>
#	include <macadam/details/numa/mc_unitnxn.h>
#	include <macadam/details/numa/mc_vander1xn.h>
#	include <macadam/details/numa/mc_var1xn.h>
#	include <macadam/details/numa/mc_zeros1x2.h>
#	include <macadam/details/numa/mc_zeros1x3.h>
#	include <macadam/details/numa/mc_zeros1xn.h>
#	include <macadam/details/numa/mc_zeros2x2.h>
#	include <macadam/details/numa/mc_zeros3x1.h>
#	include <macadam/details/numa/mc_zeros3x3.h>
#	include <macadam/details/numa/mc_zerosmx1.h>
#	include <macadam/details/numa/mc_zerosmxn.h>
#	include <macadam/details/numa/mc_zerosnxn.h>
#	include <macadam/details/numa/mc_zeye2x2.h>
#	include <macadam/details/numa/mc_zeye3x3.h>
#	include <macadam/details/numa/mc_zlogspace1xn.h>
#	include <macadam/details/numa/mc_znorm1x2.h>
#	include <macadam/details/numa/mc_znorm1x3.h>
#	include <macadam/details/numa/mc_zreig2x2.h>
#	include <macadam/details/numa/mc_zunit1x2.h>
#	include <macadam/details/numa/mc_zunit1x3.h>

#endif /* !MC_NUMA_H */

/* EOF */