# makefile for building IntegrIT(R) NatureDSP_Signal(tm) library
# Version 0.32 tuned for Cadence HiFi cores
# ---------------------------------------------------
# Target platform: gcc/xcc
#
#  make invocation:
#  make clean   - clean up all temporary files
#  make depends - refresh dependencies (makedeps.mk)
#  make all     - make depends and release
#
#  System variables XTENSA_CORE,XTENSA_SYSTEM should 
#  be properly configured for correct procedure.
#  if XTENSA_SYSTEM is undefined, g++ build procedure
#  will use cstub files and headers from local
#  directory, not from system path
#
#  Other variables:
#  LANG=C++    build under g++, otherwise xcc used
#  LANG=LLVM   build under llvm Xtensa compiler
#  SAVETEMPS=1 save temporary files
#  MEM_MODEL=1,2 select memory modeling
#  DATA_CACHE=xxx select size of data cache (applied if MEM_MODEL==2)
#-----------------------------------------------------

SRCDIRS		=	../../../../library/ \
				../../../../library/audio \
				../../../../library/complexf \
				../../../../library/diag \
				../../../../library/fft \
				../../../../library/fft/dct \
				../../../../library/fft/dct2d \
				../../../../library/fft/dctf \
				../../../../library/fft/fft \
				../../../../library/fft/fft_ie \
				../../../../library/fft/fft_ief \
				../../../../library/fft/fft_spectrum \
				../../../../library/fft/mdct \
				../../../../library/img \
				../../../../library/img/2dfft \
				../../../../library/img/imgrotate \
				../../../../library/img/imgresize \
				../../../../library/fir \
				../../../../library/fir/conv2d \
				../../../../library/fir/firblk \
				../../../../library/fir/firblkf \
				../../../../library/fir/firdec \
				../../../../library/fir/firdecf \
				../../../../library/fir/firint \
				../../../../library/fir/firintf \
				../../../../library/fir/firother \
				../../../../library/fir/firotherf \
				../../../../library/fit \
				../../../../library/fitf \
				../../../../library/iir \
				../../../../library/id \
				../../../../library/iirf \
				../../../../library/math \
				../../../../library/mathf \
				../../../../library/matop \
				../../../../library/matopf \
				../../../../library/matinv \
				../../../../library/matinv/chol32 \
				../../../../library/matinvf \
				../../../../library/matinvf/chol \
				../../../../library/tables \
				../../../../library/twiddles \
				../../../../library/vectorf \
				../../../../library/vectoref \
				../../../../library/vector 
OUTDIR		=	../../../bin
INCDIRS		=	../../../../library/include	\
			../../../../library/include_private	

TARGET				= NatureDSP_Signal
TARGET_SYMFILE			= symbols.txt
BUILD_TYPE			= LIBRARY
INTEGRIT_OLD_SUFFIXES		= 1
COMPILE_CSTUB			= 0

# check if c++ compilation required
ifeq (C++,$(LANG))
CPP	= 1
HAVE_XTENSA = 0
else
CPP	= 0
HAVE_XTENSA = 1
endif

ifeq (,$(XTENSA_SYSTEM))
ifeq (,$(CSTUB_DIR))
CSTUB_DIR :=../../../../cstub
endif
XTCORE := $(notdir $(wildcard $(CSTUB_DIR)/*-0.cpp))
XTCORE := $(subst -0.cpp,,$(XTCORE))
XTCORE := $(subst cstub-Xm_,,$(XTCORE))
else
XTCORE := $(XTENSA_CORE)
endif

ifeq (0,$(HAVE_XTENSA))
USE_PATCHES=1  # use patched cstubs for gcc, otherwise use REFERENCE variable
else
ifeq (1,$(REFERENCE))
USE_PATCHES=0
REF_SUFFIX=_reference
else
USE_PATCHES=1
REF_SUFFIX=
endif
endif

ifeq (1,$(MEM_MODEL))
MEM_SUFFIX=_mem_model
else
ifeq (2,$(MEM_MODEL))
MEM_SUFFIX=_mem_model2
else
MEM_SUFFIX=
endif
endif

ifeq (, $(DATA_CACHE))
DATA_CACHE=32768
endif

MODULES =  \
		bqriir16x16_df1_hifi4 \
		bqriir16x16_df1_nd_hifi4 \
		bqriir16x16_df2_hifi4 \
		bqriir16x16_df2_nd_hifi4 \
		bqriir32x16_df1_hifi4 \
		bqriir32x16_df1_process_hifi4 \
		bqriir32x16_df2_hifi4 \
		bqriir32x16_df2_process_hifi4 \
		bqriir32x32_df1_hifi4 \
		bqriir32x32_df1_process_hifi4 \
		bqriir32x32_df2_hifi4 \
		bqriir32x32_df2_process_hifi4 \
		bqriir32x16_df1_nd_hifi4 \
		bqriir32x16_df1_nd_process_hifi4 \
		bqriir32x16_df2_nd_hifi4 \
		bqriir32x16_df2_nd_process_hifi4 \
		bqriir32x32_df1_nd_hifi4 \
		bqriir32x32_df1_nd_process_hifi4 \
		bqriir32x32_df2_nd_hifi4 \
		bqriir32x32_df2_nd_process_hifi4 \
         scl_complex2invmag_hifi4 \
         scl_complex2mag_hifi4 \
         vec_complex2invmag_hifi4 \
         vec_complex2mag_hifi4 \
	fft_spectrumf_hifi4 \
	fft_spectrum16x32_hifi4 \
	fft_spectrum32x32_hifi4 \
	scl_add_32x16ef_hifi4 \
	scl_mul_32x16ef_hifi4 \
	scl_mac_32x16ef_hifi4 \
	vec_add_32x16ef_hifi4 \
	vec_mul_32x16ef_hifi4 \
	vec_mac_32x16ef_hifi4 \
	vec_dot_32x16ef_hifi4 \
        vec_dot_batch16x16_fast_hifi4 \
        vec_dot_batch16x16_hifi4 \
        vec_dot_batch8x16_fast_hifi4 \
        vec_dot_batch8x16_hifi4 \
        vec_dot_batch8x8_fast_hifi4 \
        vec_dot_batch8x8_hifi4 \
        vec_dot_batchf_fast_hifi4 \
        vec_dot_batchf_hifi4 \
          dct_16x16_cffts_hifi4 \
          dct_16x16_hifi4 \
          dct_32x16_cffts_hifi4 \
          dct_32x16_hifi4 \
          dct_32x32_cffts_hifi4 \
          dct_32x32_hifi4 \
          dct4_32x16_hifi4 \
          dct4_32x32_hifi4 \
          dct2d_8x16_hifi4 \
          idct2d_16x8_hifi4 \
          dctf_hifi4 \
          fft_cplx_16x16_hifi4 \
          fft_cplx_16x16_stages_scl2_r2_r3_r5_hifi4 \
          fft_cplx_16x16_stages_scl2_r4_r8_hifi4 \
          fft_cplx_16x16_stages_scl3_r2_r3_r5_hifi4 \
          fft_cplx_16x16_stages_scl3_r4_r8_hifi4 \
          fft_cplx_32x16_hifi4 \
          fft_cplx_32x16_stages_scl2_r2_r3_r5_hifi4 \
          fft_cplx_32x16_stages_scl2_r4_r8_hifi4 \
          fft_cplx_32x16_stages_scl3_r2_r3_r5_hifi4 \
          fft_cplx_32x16_stages_scl3_r4_r8_hifi4 \
          fft_cplx_32x32_hifi4 \
          fft_cplx_stages_S2_radix2_3_5_32x32_hifi4 \
          fft_cplx_stages_S2_radix4_8_32x32_hifi4 \
          fft_cplx_stages_S3_radix2_3_5_32x32_hifi4 \
          fft_cplx_stages_S3_radix4_8_32x32_hifi4 \
          fft_real_16x16_hifi4 \
          fft_real_32x16_hifi4 \
          fft_real_32x32_hifi4 \
          ifft_cplx_16x16_hifi4 \
          ifft_cplx_32x16_hifi4 \
          ifft_cplx_32x32_hifi4 \
          ifft_real_16x16_hifi4 \
          ifft_real_32x16_hifi4 \
          ifft_real_32x32_hifi4 \
          fft_cplx16x16_ie_hifi4 \
          fft_cplx32x16_ie_hifi4 \
          fft_cplx32x32_ie_hifi4 \
          fft_real16x16_ie_hifi4 \
          fft_real32x16_ie_hifi4 \
          fft_real32x32_ie_hifi4 \
          fft_stage_inner_DFT4_16x16_ie_hifi4 \
          ifft_cplx16x16_ie_hifi4 \
          ifft_cplx32x16_ie_hifi4 \
          ifft_cplx32x32_ie_hifi4 \
          ifft_real16x16_ie_hifi4 \
          ifft_real32x16_ie_hifi4 \
          ifft_real32x32_ie_hifi4 \
          fft_cplxf_ie_hifi4 \
          fft_realf_ie_hifi4 \
          ifft_cplxf_ie_hifi4 \
          ifft_realf_ie_hifi4 \
          stereo_fft_cplx16x16_ie_hifi4 \
          stereo_fft_cplx32x16_ie_hifi4 \
          stereo_fft_cplx32x32_ie_hifi4 \
          stereo_fft_cplxf_ie_hifi4 \
          stereo_ifft_cplx16x16_ie_hifi4 \
          stereo_ifft_cplx32x16_ie_hifi4 \
          stereo_ifft_cplx32x32_ie_hifi4 \
          stereo_ifft_cplxf_ie_hifi4 \
          stereo_fft_cplx32x16_ie_inner_hifi4 \
          stereo_fft_cplx32x32_ie_inner_hifi4 \
          imdct_32x16_hifi4 \
          imdct_32x32_hifi4 \
          mdct_32x16_hifi4 \
          mdct_32x32_hifi4 \
          conv2d_11x7_16x16_hifi4 \
          conv2d_11x7_8x16_hifi4 \
          conv2d_11x7_8x8_hifi4 \
          conv2d_3x3_16x16_hifi4 \
          conv2d_3x3_8x16_hifi4 \
          conv2d_3x3_8x8_hifi4 \
          conv2d_5x5_16x16_hifi4 \
          conv2d_5x5_8x16_hifi4 \
          conv2d_5x5_8x8_hifi4 \
          conv2d_11x7f_hifi4 \
          conv2d_3x3f_hifi4 \
          conv2d_5x5f_hifi4 \
       		 stereo_bkfir16x16_hifi4 \
		 stereo_bkfir16x16_process_hifi4 \
		 stereo_bkfir32x32_hifi4 \
		 stereo_bkfir32x32_process_hifi4 \
		 stereo_bqriir16x16_df1_hifi4 \
		 stereo_bqriir16x16_df1_process_hifi4 \
		 stereo_bqriir32x16_df1_hifi4 \
		 stereo_bqriir32x16_df1_process_hifi4 \
		 stereo_bqriir32x32_df1_hifi4 \
		 stereo_bqriir32x32_df1_process_hifi4 \
		 stereo_bqriirf_df1_hifi4 \
		 stereo_bqriirf_df1_process_hifi4 \
		 stereo_bqriir16x16_df1_nd_hifi4 \
		 stereo_bqriir16x16_df1_nd_process_hifi4 \
		 stereo_bqriir32x16_df1_nd_hifi4 \
		 stereo_bqriir32x16_df1_nd_process_hifi4 \
		 stereo_bqriir32x32_df1_nd_hifi4 \
		 stereo_bqriir32x32_df1_nd_process_hifi4 \
		 stereo_bqriirf_df1_nd_hifi4 \
		 stereo_bqriirf_df1_nd_process_hifi4 \
		 stereo_bkfirf_hifi4 \
		 stereo_bkfirf_process_hifi4 \
       		 bkfir24x24p_hifi4 \
       		 bkfir24x24p_process_hifi4 \
       		 bkfir16x16_hifi4 \
		 bkfir32x16_hifi4 \
		 bkfir32x16_process_hifi4 \
		 bkfir32x32ep_hifi4 \
		 bkfir32x32ep_process_hifi4 \
		 bkfir32x32_hifi4 \
		 bkfir32x32_process_hifi4 \
		 bkfira16x16_hifi4 \
		 bkfira32x16_hifi4 \
		 bkfira32x16_process_hifi4 \
		 bkfira32x32ep_hifi4 \
		 bkfira32x32ep_process_hifi4 \
		 bkfira32x32_hifi4 \
		 bkfira32x32_process_hifi4 \
		 cxfir16x16_hifi4 \
		 cxfir32x16_hifi4 \
		 cxfir32x16_process_hifi4 \
		 cxfir32x32ep_hifi4 \
		 cxfir32x32ep_process_hifi4 \
		 cxfir32x32_hifi4 \
		 cxfir32x32_process_hifi4 \
		 bkfiraf_hifi4 \
		 bkfiraf_process_hifi4 \
		 bkfirf_hifi4 \
		 bkfirf_process_hifi4 \
		 cxfirf_hifi4 \
		 cxfirf_process_hifi4 \
		 firdec16x16_D2_hifi4 \
		 firdec16x16_D3_hifi4 \
		 firdec16x16_D4_hifi4 \
		 firdec16x16_DX_hifi4 \
		 firdec16x16_hifi4 \
		 firdec32x16_D2_hifi4 \
		 firdec32x16_D3_hifi4 \
		 firdec32x16_D4_hifi4 \
		 firdec32x16_DX_hifi4 \
		 firdec32x16_hifi4 \
		 firdec32x32ep_D2_hifi4 \
		 firdec32x32ep_D3_hifi4 \
		 firdec32x32ep_D4_hifi4 \
		 firdec32x32ep_DX_hifi4 \
		 firdec32x32ep_hifi4 \
		 firdec32x32_D2_hifi4 \
		 firdec32x32_D3_hifi4 \
		 firdec32x32_D4_hifi4 \
		 firdec32x32_DX_hifi4 \
		 firdec32x32_hifi4 \
		 firdecf_hifi4 \
		 fir_decimaf_2x_hifi4 \
		 fir_decimaf_3x_hifi4 \
		 fir_decimaf_4x_hifi4 \
		 fir_decimaf_Dx_hifi4 \
		 firinterp16x16_D2_hifi4 \
		 firinterp16x16_D3_hifi4 \
		 firinterp16x16_D4_hifi4 \
		 firinterp16x16_DX_hifi4 \
		 firinterp16x16_hifi4 \
		 firinterp32x16_D2_hifi4 \
		 firinterp32x16_D3_hifi4 \
		 firinterp32x16_D4_hifi4 \
		 firinterp32x16_DX_hifi4 \
		 firinterp32x16_hifi4 \
		 firinterp32x32ep_D2_hifi4 \
		 firinterp32x32ep_D3_hifi4 \
		 firinterp32x32ep_D4_hifi4 \
		 firinterp32x32ep_DX_hifi4 \
		 firinterp32x32ep_hifi4 \
		 firinterp32x32_D2_hifi4 \
		 firinterp32x32_D3_hifi4 \
		 firinterp32x32_D4_hifi4 \
		 firinterp32x32_DX_hifi4 \
		 firinterp32x32_hifi4 \
		 firinterpf_hifi4 \
		 fir_interpf_2x_hifi4 \
		 fir_interpf_3x_hifi4 \
		 fir_interpf_4x_hifi4 \
		 fir_interpf_Dx_hifi4 \
		 cxfir_convol32x16_hifi4 \
		 cxfir_convola32x16_hifi4 \
		 fir_acorr16x16_hifi4 \
		 fir_acorr32x32ep_hifi4 \
		 fir_acorr32x32_hifi4 \
		 fir_acorra16x16_hifi4 \
		 fir_acorra32x32ep_hifi4 \
		 fir_acorra32x32_hifi4 \
		 fir_blms16x16_hifi4 \
		 fir_blms16x32_hifi4 \
		 fir_blms32x32ep_hifi4 \
		 fir_blms32x32_hifi4 \
		 cxfir_blms32x32_hifi4 \
		 fir_convol16x16_hifi4 \
		 fir_convol32x16_hifi4 \
		 fir_convol32x32ep_hifi4 \
		 fir_convol32x32_hifi4 \
		 fir_convola16x16_hifi4 \
		 fir_convola32x16_hifi4 \
		 fir_convola32x32ep_hifi4 \
		 fir_convola32x32_hifi4 \
		 fir_lacorra16x16_hifi4 \
		 fir_lacorra32x32_hifi4 \
		 fir_lconvola16x16_hifi4 \
		 fir_lconvola32x32_hifi4 \
		 fir_lxcorra16x16_hifi4 \
		 fir_lxcorra32x32_hifi4 \
		 fir_xcorr16x16_hifi4 \
		 fir_xcorr32x16_hifi4 \
		 fir_xcorr32x32ep_hifi4 \
		 fir_xcorr32x32_hifi4 \
		 cxfir_xcorr32x32_hifi4 \
		 cxfir_xcorra32x32_hifi4 \
		 fir_xcorra16x16_hifi4 \
		 fir_xcorra32x16_hifi4 \
		 fir_xcorra32x32ep_hifi4 \
		 fir_xcorra32x32_hifi4 \
		 raw_corr16x16_hifi4 \
		 raw_corr32x16_hifi4 \
		 raw_corr32x32_hifi4 \
		 raw_corr32x32ep_hifi4 \
		 raw_lxcorr16x16_hifi4 \
		 raw_lxcorr32x32_hifi4 \
		 cxfir_xcorraf_hifi4 \
		 cxfir_xcorrf_hifi4 \
		 fir_acorraf_hifi4 \
		 fir_acorrf_hifi4 \
		 fir_blmsf_hifi4 \
		 cxfir_blmsf_hifi4 \
		 fir_convolaf_hifi4 \
		 fir_convolf_hifi4 \
		 fir_xcorraf_hifi4 \
		 fir_xcorrf_hifi4 \
		 raw_corrf_hifi4 \
		 vec_poly4_32x32_hifi4 \
		 vec_poly8_32x32_hifi4 \
		 vec_poly4f_hifi4 \
		 vec_poly8f_hifi4 \
		latr16x16_hifi4 \
		latr32x16_1_proc_hifi4 \
		latr32x16_2_proc_hifi4 \
		latr32x16_3_proc_hifi4 \
		latr32x16_4_proc_hifi4 \
		latr32x16_5_proc_hifi4 \
		latr32x16_6_proc_hifi4 \
		latr32x16_7_proc_hifi4 \
		latr32x16_8_proc_hifi4 \
		latr32x16_hifi4 \
		latr32x16_X_proc_hifi4 \
		latr32x32_1_proc_hifi4 \
		latr32x32_2_proc_hifi4 \
		latr32x32_3_proc_hifi4 \
		latr32x32_4_proc_hifi4 \
		latr32x32_5_proc_hifi4 \
		latr32x32_6_proc_hifi4 \
		latr32x32_7_proc_hifi4 \
		latr32x32_8_proc_hifi4 \
		latr32x32_hifi4 \
		latr32x32_X_proc_hifi4 \
		bqciirf_df1_hifi4 \
		bqciirf_df1_process_hifi4 \
		bqriirf_df1_hifi4 \
		bqriirf_df1_process_hifi4 \
		bqriirf_df2t_hifi4 \
		bqriirf_df2t_process_hifi4 \
		bqriirf_df2_hifi4 \
		bqriirf_df2_process_hifi4 \
		bqciirf_df1_nd_hifi4 \
		bqciirf_df1_nd_process_hifi4 \
		bqriirf_df1_nd_hifi4 \
		bqriirf_df1_nd_process_hifi4 \
		bqriirf_df2t_nd_hifi4 \
		bqriirf_df2t_nd_process_hifi4 \
		bqriirf_df2_nd_hifi4 \
		bqriirf_df2_nd_process_hifi4 \
		latrf1_hifi4 \
		latrf2_hifi4 \
		latrf3_hifi4 \
		latrf4_hifi4 \
		latrf5_hifi4 \
		latrf6_hifi4 \
		latrf7_hifi4 \
		latrf8_hifi4 \
		latrfX_hifi4 \
		latrf_hifi4 \
		scl_alog10_32x32_hifi4 \
		scl_alog2_32x32_hifi4 \
		scl_alogn_32x32_hifi4 \
		scl_atan_32x32_hifi4 \
		scl_cosine_32x32_hifi4 \
		scl_divide16x16_hifi4 \
		scl_divide32x32_hifi4 \
		scl_divide64x32_hifi4 \
		scl_log10_32x32_hifi4 \
		scl_log2_32x32_hifi4 \
		scl_logn_32x32_hifi4 \
		scl_recip16x16_hifi4 \
		scl_recip32x32_hifi4 \
		scl_rsqrt16x16_hifi4 \
		scl_rsqrt32x32_hifi4 \
		scl_sigmoid32x32_hifi4 \
		scl_sine_32x32_hifi4 \
		scl_sqrt16x16_hifi4 \
		scl_sqrt64x32_hifi4 \
		scl_sqrt_32x32_hifi4 \
		scl_tanh32x32_hifi4 \
		scl_tan_32x32_hifi4 \
		vec_alog10_32x32_hifi4 \
		vec_alog2_32x32_hifi4 \
		vec_alogn_32x32_hifi4 \
		vec_atan_32x32_hifi4 \
		vec_cosine_32x32_fast_hifi4 \
		vec_cosine_32x32_hifi4 \
		vec_divide16x16_fast_hifi4 \
		vec_divide16x16_hifi4 \
		vec_divide32x32_fast_hifi4 \
		vec_divide32x32_hifi4 \
		vec_divide64x32i_hifi4 \
		vec_log10_32x32_hifi4 \
		vec_log2_32x32_hifi4 \
		vec_logn_32x32_hifi4 \
		vec_recip16x16_hifi4 \
		vec_recip32x32_hifi4 \
		vec_rsqrt16x16_hifi4 \
		vec_rsqrt32x32_hifi4 \
		vec_sigmoid32x32_hifi4 \
		vec_sine_32x32_fast_hifi4 \
		vec_sine_32x32_hifi4 \
		vec_softmax32x32_hifi4 \
		vec_sqrt16x16_hifi4 \
		vec_sqrt64x32_hifi4 \
		vec_sqrt_32x32_fast_hifi4 \
		vec_sqrt_32x32_hifi4 \
		vec_tanh32x32_hifi4 \
		vec_tan_32x32_hifi4 \
		scl_antilog10f_hifi4 \
		scl_antilog2f_hifi4 \
		scl_antilognf_hifi4 \
		scl_atan2f_hifi4 \
		scl_atanf_hifi4 \
		scl_cosinef_hifi4 \
		scl_float2int_hifi4 \
		scl_int2float_hifi4 \
		scl_log10f_hifi4 \
		scl_log2f_hifi4 \
		scl_lognf_hifi4 \
		scl_sinef_hifi4 \
		scl_tanf_hifi4 \
		scl_tanhf_hifi4 \
		scl_reluf_hifi4 \
		scl_relu32x32_hifi4 \
		scl_recip64x64_hifi4 \
		scl_divide64x64_hifi4 \
		scl_sqrt32x16_hifi4 \
		scl_sigmoidf_hifi4 \
		vec_alog10f_hifi4 \
		vec_alog2f_hifi4 \
		vec_alognf_hifi4 \
		vec_atan2f_hifi4 \
		vec_atanf_hifi4 \
		vec_cosinef_hifi4 \
		vec_float2int_hifi4 \
		vec_int2float_hifi4 \
		vec_log10f_hifi4 \
		vec_log2f_hifi4 \
		vec_lognf_hifi4 \
		vec_sinef_hifi4 \
		vec_tanf_hifi4 \
		vec_tanhf_hifi4 \
		vec_reluf_hifi4 \
		vec_relu32x32_hifi4 \
		vec_softmaxf_hifi4 \
		vec_recip64x64_hifi4 \
		vec_divide64x64_hifi4 \
		vec_sqrt32x16_hifi4 \
		vec_sigmoidf_hifi4 \
		vec_pow_32x32_hifi4 \
		scl_pow_32x32_hifi4 \
		vec_powf_hifi4 \
		scl_powf_hifi4 \
		mtx_inv2x2f_hifi4 \
		mtx_inv3x3f_hifi4 \
		mtx_inv4x4f_hifi4 \
		mtx_inv6x6f_hifi4 \
		mtx_inv8x8f_hifi4 \
		mtx_inv10x10f_hifi4 \
		cmtx_gjelim10x10_32x32_hifi4 \
		cmtx_gjelim2x2_32x32_hifi4 \
		cmtx_gjelim3x3_32x32_hifi4 \
		cmtx_gjelim4x4_32x32_hifi4 \
		cmtx_gjelim6x6_32x32_hifi4 \
		cmtx_gjelim8x8_32x32_hifi4 \
		cmtx_inv10x10_32x32_hifi4 \
		cmtx_inv2x2_32x32_hifi4 \
		cmtx_inv3x3_32x32_hifi4 \
		cmtx_inv4x4_32x32_hifi4 \
		cmtx_inv6x6_32x32_hifi4 \
		cmtx_inv8x8_32x32_hifi4 \
		mtx_gjelim10x10_32x32_hifi4 \
		mtx_gjelim2x2_32x32_hifi4 \
		mtx_gjelim3x3_32x32_hifi4 \
		mtx_gjelim4x4_32x32_hifi4 \
		mtx_gjelim6x6_32x32_hifi4 \
		mtx_gjelim8x8_32x32_hifi4 \
		mtx_inv10x10_32x32_hifi4 \
		mtx_inv2x2_32x32_hifi4 \
		mtx_inv3x3_32x32_hifi4 \
		mtx_inv4x4_32x32_hifi4 \
		mtx_inv6x6_32x32_hifi4 \
		mtx_inv8x8_32x32_hifi4 \
		mtx_mpy16x16_fast_hifi4 \
		mtx_mpy16x16_hifi4 \
		mtx_mpy32x32_fast_hifi4 \
		mtx_mpy32x32_hifi4 \
		mtx_mpy8x16_fast_hifi4 \
		mtx_mpy8x16_hifi4 \
		mtx_mpy8x8_fast_hifi4 \
		mtx_mpy8x8_hifi4 \
		mtx_mpyt16x16_fast_hifi4 \
		mtx_mpyt16x16_hifi4 \
		mtx_mpyt32x32_fast_hifi4 \
		mtx_mpyt32x32_hifi4 \
		mtx_mpyt8x16_fast_hifi4 \
		mtx_mpyt8x16_hifi4 \
		mtx_mpyt8x8_fast_hifi4 \
		mtx_mpyt8x8_hifi4 \
		mtx_transpose16x16_fast_hifi4 \
		mtx_transpose16x16_hifi4 \
		mtx_transpose32x32_fast_hifi4 \
		mtx_transpose32x32_hifi4 \
		mtx_transpose8x8_fast_hifi4 \
		mtx_transpose8x8_hifi4 \
		mtx_transposef_fast_hifi4 \
		mtx_transposef_hifi4 \
		mtx_vecmpy16x16_fast_hifi4 \
		mtx_vecmpy16x16_hifi4 \
		mtx_vecmpy32x32_fast_hifi4 \
		mtx_vecmpy32x32_hifi4 \
		mtx_vecmpy8x16_fast_hifi4 \
		mtx_vecmpy8x16_hifi4 \
		mtx_vecmpy8x8_fast_hifi4 \
		mtx_vecmpy8x8_hifi4 \
		mtx_mpyf_fast_hifi4 \
		mtx_mpyf_hifi4 \
		mtx_mpytf_fast_hifi4 \
		mtx_mpytf_hifi4 \
		mtx_vecmpyf_fast_hifi4 \
		mtx_vecmpyf_hifi4 \
		scl_bexp16_hifi4 \
		scl_bexp32_hifi4 \
		vec_add16x16_fast_hifi4 \
		vec_add16x16_hifi4 \
		vec_add32x32_fast_hifi4 \
		vec_add32x32_hifi4 \
		vec_bexp16_fast_hifi4 \
		vec_bexp16_hifi4 \
		vec_bexp32_fast_hifi4 \
		vec_bexp32_hifi4 \
		vec_dot16x16_fast_hifi4 \
		vec_dot16x16_hifi4 \
		vec_dot32x16_fast_hifi4 \
		vec_dot32x16_hifi4 \
		vec_dot32x32_fast_hifi4 \
		vec_dot32x32_hifi4 \
		vec_dot64x32_fast_hifi4 \
		vec_dot64x32_hifi4 \
		vec_dot64x64i_fast_hifi4 \
		vec_dot64x64i_hifi4 \
		vec_dot64x64_fast_hifi4 \
		vec_dot64x64_hifi4 \
		vec_max_16x16_fast_hifi4 \
		vec_max_16x16_hifi4 \
		vec_max_32x32_fast_hifi4 \
		vec_max_32x32_hifi4 \
		vec_min_16x16_fast_hifi4 \
		vec_min_16x16_hifi4 \
		vec_min_32x32_fast_hifi4 \
		vec_min_32x32_hifi4 \
		vec_power16x16_fast_hifi4 \
		vec_power16x16_hifi4 \
		vec_power32x32_fast_hifi4 \
		vec_power32x32_hifi4 \
		vec_scale16x16_fast_hifi4 \
		vec_scale16x16_hifi4 \
		vec_scale32x32_fast_hifi4 \
		vec_scale32x32_hifi4 \
		vec_shift16x16_fast_hifi4 \
		vec_shift16x16_hifi4 \
		vec_shift32x32_fast_hifi4 \
		vec_shift32x32_hifi4 \
		scl_bexpf_hifi4 \
		vec_addf_hifi4 \
		vec_bexpf_hifi4 \
		vec_dotf_hifi4 \
		vec_maxf_hifi4 \
		vec_minf_hifi4 \
		vec_powerf_hifi4 \
		vec_scalef_hifi4 \
		vec_scale_sf_hifi4 \
		vec_shiftf_hifi4 \
		alog10f_tbl \
		alog2f_tbl \
		atanf_tbl \
		expf_tbl \
		inff_tbl \
		inv2pif_tbl \
		log10f_tbl \
		log2f_tbl \
		lognf_tbl \
		nan_tbl \
		nanf_tbl \
		pif_tbl \
		polyrsqrtq23_tbl \
		scl_atan_table \
		scl_atan_table16 \
         scl_sine_table16 \
         scl_sine_table32 \
         scl_sqrt_table \
         sinf_tbl \
         sqrt2f_tbl \
         tanf_tbl \
	tanhf_tbl \
	pow2f_tbl \
         vec_alog_table \
         vec_log_table \
         vec_recip_table \
         vec_pow_32x32_table \
	logmel32x32_hifi4 \
	logmelf_hifi4 \
	mfcc32x32_hifi4 \
	mfccf_hifi4 \
          dct_16_32 \
          dct_16_64 \
          dct_32_32 \
          dct_32_64 \
          dct4_16_128 \
          dct4_16_256 \
          dct4_16_32 \
          dct4_16_512 \
          dct4_16_64 \
          dct4_32_128 \
          dct4_32_256 \
          dct4_32_32 \
          dct4_32_512 \
          dct4_32_64 \
          dctf_32_twd \
          dctf_64_twd \
          fft_cplx_twd16_1024_tbl \
          fft_cplx_twd16_128_tbl \
          fft_cplx_twd16_144_tbl \
          fft_cplx_twd16_16_tbl \
          fft_cplx_twd16_160_tbl \
          fft_cplx_twd16_176_tbl \
          fft_cplx_twd16_192_tbl \
          fft_cplx_twd16_2048_tbl \
          fft_cplx_twd16_240_tbl \
          fft_cplx_twd16_256_tbl \
          fft_cplx_twd16_288_tbl \
          fft_cplx_twd16_32_tbl \
          fft_cplx_twd16_320_tbl \
          fft_cplx_twd16_352_tbl \
          fft_cplx_twd16_384_tbl \
          fft_cplx_twd16_4096_tbl \
          fft_cplx_twd16_48_tbl \
          fft_cplx_twd16_480_tbl \
          fft_cplx_twd16_512_tbl \
          fft_cplx_twd16_576_tbl \
          fft_cplx_twd16_64_tbl \
          fft_cplx_twd16_640_tbl \
          fft_cplx_twd16_96_tbl \
          fft_cplx_twd32_100_tbl \
          fft_cplx_twd32_1024_tbl \
          fft_cplx_twd32_108_tbl \
          fft_cplx_twd32_12_tbl \
          fft_cplx_twd32_120_tbl \
          fft_cplx_twd32_128_tbl \
          fft_cplx_twd32_144_tbl \
          fft_cplx_twd32_16_tbl \
          fft_cplx_twd32_160_tbl \
          fft_cplx_twd32_180_tbl \
          fft_cplx_twd32_192_tbl \
          fft_cplx_twd32_200_tbl \
          fft_cplx_twd32_2048_tbl \
          fft_cplx_twd32_216_tbl \
          fft_cplx_twd32_24_tbl \
          fft_cplx_twd32_240_tbl \
          fft_cplx_twd32_256_tbl \
          fft_cplx_twd32_288_tbl \
          fft_cplx_twd32_300_tbl \
          fft_cplx_twd32_32_tbl \
          fft_cplx_twd32_320_tbl \
          fft_cplx_twd32_324_tbl \
          fft_cplx_twd32_36_tbl \
          fft_cplx_twd32_360_tbl \
          fft_cplx_twd32_384_tbl \
          fft_cplx_twd32_400_tbl \
          fft_cplx_twd32_4096_tbl \
          fft_cplx_twd32_432_tbl \
          fft_cplx_twd32_48_tbl \
          fft_cplx_twd32_480_tbl \
          fft_cplx_twd32_512_tbl \
          fft_cplx_twd32_540_tbl \
          fft_cplx_twd32_576_tbl \
          fft_cplx_twd32_60_tbl \
          fft_cplx_twd32_600_tbl \
          fft_cplx_twd32_64_tbl \
          fft_cplx_twd32_72_tbl \
          fft_cplx_twd32_768_tbl \
          fft_cplx_twd32_80_tbl \
          fft_cplx_twd32_96_tbl \
          fft_cplx_twd32_960_tbl \
          fft_cplx_twd32x16_160_tbl \
          fft_cplx_twd32x16_192_tbl \
          fft_cplx_twd32x16_240_tbl \
          fft_cplx_twd32x16_320_tbl \
          fft_cplx_twd32x16_384_tbl \
          fft_cplx_twd32x16_480_tbl \
          fft_real_twd16_1024_tbl \
          fft_real_twd16_128_tbl \
          fft_real_twd16_144_tbl \
          fft_real_twd16_160_tbl \
          fft_real_twd16_176_tbl \
          fft_real_twd16_192_tbl \
          fft_real_twd16_2048_tbl \
          fft_real_twd16_240_tbl \
          fft_real_twd16_256_tbl \
          fft_real_twd16_288_tbl \
          fft_real_twd16_32_tbl \
          fft_real_twd16_320_tbl \
          fft_real_twd16_352_tbl \
          fft_real_twd16_384_tbl \
          fft_real_twd16_4096_tbl \
          fft_real_twd16_480_tbl \
          fft_real_twd16_512_tbl \
          fft_real_twd16_576_tbl \
          fft_real_twd16_64_tbl \
          fft_real_twd16_640_tbl \
          fft_real_twd16_8192_tbl \
          fft_real_twd16_96_tbl \
          fft_real_twd32_1024_tbl \
          fft_real_twd32_108_tbl \
          fft_real_twd32_1152_tbl \
          fft_real_twd32_12_tbl \
          fft_real_twd32_120_tbl \
          fft_real_twd32_128_tbl \
          fft_real_twd32_144_tbl \
          fft_real_twd32_1440_tbl \
          fft_real_twd32_1536_tbl \
          fft_real_twd32_160_tbl \
          fft_real_twd32_180_tbl \
          fft_real_twd32_192_tbl \
          fft_real_twd32_1920_tbl \
          fft_real_twd32_2048_tbl \
          fft_real_twd32_216_tbl \
          fft_real_twd32_24_tbl \
          fft_real_twd32_240_tbl \
          fft_real_twd32_256_tbl \
          fft_real_twd32_288_tbl \
          fft_real_twd32_30_tbl \
          fft_real_twd32_300_tbl \
          fft_real_twd32_32_tbl \
          fft_real_twd32_320_tbl \
          fft_real_twd32_324_tbl \
          fft_real_twd32_36_tbl \
          fft_real_twd32_360_tbl \
          fft_real_twd32_384_tbl \
          fft_real_twd32_4096_tbl \
          fft_real_twd32_432_tbl \
          fft_real_twd32_48_tbl \
          fft_real_twd32_480_tbl \
          fft_real_twd32_512_tbl \
          fft_real_twd32_540_tbl \
          fft_real_twd32_576_tbl \
          fft_real_twd32_60_tbl \
          fft_real_twd32_64_tbl \
          fft_real_twd32_72_tbl \
          fft_real_twd32_720_tbl \
          fft_real_twd32_768_tbl \
          fft_real_twd32_8192_tbl \
          fft_real_twd32_90_tbl \
          fft_real_twd32_96_tbl \
          fft_real_twd32_960_tbl \
          fft_real_twd32x16_160_tbl \
          fft_real_twd32x16_192_tbl \
          fft_real_twd32x16_240_tbl \
          fft_real_twd32x16_320_tbl \
          fft_real_twd32x16_384_tbl \
          fft_real_twd32x16_480_tbl \
         cubic_kernel_hifi4 \
         divide_q31_hifi4 \
         img_getCoef_dn_cubic_hifi4 \
         img_getCoef_dn_hifi4 \
         img_getCoef_up_cubic_hifi4 \
         img_getCoef_up_hifi4 \
         imghist_gu8_hifi4 \
         imghist_gs8_hifi4 \
         imghist_gs8_hifi4 \
         img_resize_hifi4 \
         imgconvert_rgbyuv_hifi4 \
         imgconvert_rgbyuv16_hifi4 \
         imgconvert_yuvrgb_hifi4 \
         imgconvert_yuvrgb16_hifi4 \
         imgcopy_api_gu8_temp_hifi4 \
         imgcopy_api_gs8_temp_hifi4 \
         imgcopy_api_gs16_temp_hifi4 \
         imgcopy_api_gs16fast_temp_hifi4 \
         imgcopy_api_gu8fast_temp_hifi4 \
         imgcopy_api_gs8fast_temp_hifi4 \
         imgcopy_api_temp_gu8_hifi4 \
         imgcopy_api_temp_gs8_hifi4 \
         imgcopy_api_temp_gs16_hifi4 \
         imgcopy_api_temp_gs16fast_hifi4 \
         imgcopy_api_temp_gu8fast_hifi4 \
         imgcopy_api_temp_gs8fast_hifi4 \
         imgcopyHorz_gu8_hifi4 \
         imgcopyHorz_gs8_hifi4 \
         imgcopyHorz_gs16_hifi \
         imgcopyLeftRightBorders_gu8_hifi4 \
         imgcopyLeftRightBorders_gs8_hifi4 \
         imgcopyLeftRightBorders_gs16_hifi4 \
         imgcopyUpBtmBorders_gu8_hifi4 \
         imgcopyUpBtmBorders_gs8_hifi4 \
         imgcopyUpBtmBorders_gs16_hifi4 \
         imgcopyVert_gu8_hifi4 \
         imgcopyVert_gs8_hifi4 \
         imgcopyVert_gs16_hifi4 \
         imgdeinterleave_hifi4 \
         imgdeinterleave16_hifi4 \
         imgfastconvert_rgbyuv_hifi4 \
         imgfastconvert_rgbyuv16_hifi4 \
         imgfastconvert_yuvrgb_hifi4 \
         imgfastconvert_yuvrgb16_hifi4 \
         imgfastdeinterleave_hifi4 \
         imgfastdeinterleave16_hifi4 \
         imgfasthist_gu8_hifi4 \
         imgfasthist_gs8_hifi4 \
         imgfasthist_gs16_hifi4 \
         imgfastinterleave_hifi4 \
         imgfastinterleave16_hifi4 \
         imgfastnorm_gu8_hifi4 \
         imgfastnorm_gu8_nonlinear_hifi4 \
         imgfastnorm_gs8_hifi4 \
         imgfastnorm_gs8_nonlinear_hifi4 \
         imgfastnorm_gs16_hifi4 \
         imgfastnorm_gs16_nonlinear_hifi4 \
         imgfastpad_gu8_hifi4 \
         imgfastpad_gs8_hifi4 \
         imgfastpad_gs16_hifi4 \
         imgfastresize_gu8_hifi4 \
         imgfastresize_gs8_hifi4 \
         imgfastresize_gs16_hifi4 \
         imghist_gs16_hifi4 \
         imginterleave_hifi4 \
         imginterleave16_hifi4 \
         imgnorm_gu8_hifi4 \
         imgnorm_gu8_nonlinear_hifi4 \
         imgnorm_gs8_hifi4 \
         imgnorm_gs8_nonlinear_hifi4 \
         imgnorm_gs16_hifi4 \
         imgnorm_gs16_nonlinear_hifi4 \
         imgpad_gu8_hifi4 \
         imgpad_gs8_hifi4 \
         imgpad_gs16_hifi4 \
         imgpad_process_hifi4 \
         imgresize_gu8_hifi4 \
         imgresize_gs8_hifi4 \
         imgresize_gs16_hifi4 \
         imgresizer_bicubic_hifi4 \
         imgresizer_bilinear_hifi4 \
         imgresizer_dn2xh_cubic_hifi4 \
         imgresizer_dn2xh_hifi4 \
         imgresizer_dn2xv_cubic_hifi4 \
         imgresizer_dn2xv_hifi4 \
         imgresizer_dnxh_cubic_hifi4 \
         imgresizer_dnxh_hifi4 \
         imgresizer_dnxv_cubic_hifi4 \
         imgresizer_dnxv_hifi4 \
         imgresizer_nearest_hifi4 \
         imgresizer_nh_hifi4 \
         imgresizer_nv_hifi4 \
         imgresizer_up2xh_cubic_hifi4 \
         imgresizer_up2xh_hifi4 \
         imgresizer_up2xv_cubic_hifi4 \
         imgresizer_up2xv_hifi4 \
         imgresizer_upxh_cubic_hifi4 \
         imgresizer_upxh_hifi4 \
         imgresizer_upxv_cubic_hifi4 \
         imgresizer_upxv_hifi4 \
         imgsubcopy_gu8_hifi4 \
         imgsubcopy_gs8_hifi4 \
         imgsubcopy_gs16_hifi4 \
         imgfft_common_hifi4 \
         imgfft_gu8_hifi4 \
         imgfft_gs8_hifi4 \
         imgfft_gs16_hifi4 \
         imgifft_common_hifi4 \
         imgifft_gu8_hifi4 \
         imgifft_gs8_hifi4 \
         imgifft_gs16_hifi4 \
         img_gu8_crop_hifi4 \
         img_gs8_crop_hifi4 \
         img_gs_hshear_hifi4 \
         img_gu8_rot0_16_hifi4 \
         img_gu8_rot0_hifi4 \
         img_gu8_rot180_16_hifi4 \
         img_gu8_rot180_hifi4 \
         img_gu8_rot270_16_hifi4 \
         img_gu8_rot270_hifi4 \
         img_gu8_rot90_16_hifi4 \
         img_gu8_rot90_hifi4 \
         img_gu8_rotq_hifi4 \
         img_gu8_rotq16_hifi4 \
         img_gs8_rot0_16_hifi4 \
         img_gs8_rot180_16_hifi4 \
         img_gs8_rot270_16_hifi4 \
         img_gs8_rot90_16_hifi4 \
         img_gs8_rotq_hifi4 \
         img_gs8_rotq16_hifi4 \
         img_gs_vshear_hifi4 \
         img_gs16_crop_hifi4 \
         img_gs16_rot0_hifi4 \
         img_gs16_rot0_16_hifi4 \
         img_gs16_rot180_hifi4 \
         img_gs16_rot180_16_hifi4 \
         img_gs16_rot270_hifi4 \
         img_gs16_rot270_16_hifi4 \
         img_gs16_rot90_hifi4 \
         img_gs16_rot90_16_hifi4 \
         img_gs16_rotq_hifi4 \
         img_gs16_rotq16_hifi4 \
         imgfast_gu8_crop_hifi4 \
         imgfast_gu8_rot0_16_hifi4 \
         imgfast_gu8_rot0_hifi4 \
         imgfast_gu8_rot180_16_hifi4 \
         imgfast_gu8_rot180_hifi4 \
         imgfast_gu8_rot270_16_hifi4 \
         imgfast_gu8_rot270_hifi4 \
         imgfast_gu8_rot90_16_hifi4 \
         imgfast_gu8_rot90_hifi4 \
         imgfast_gu8_rotq_hifi4 \
         imgfast_gu8_rotq16_hifi4 \
         imgfast_gs8_crop_hifi4 \
         imgfast_gs8_rot0_16_hifi4 \
         imgfast_gs8_rot180_16_hifi4 \
         imgfast_gs8_rot270_16_hifi4 \
         imgfast_gs8_rot90_16_hifi4 \
         imgfast_gs8_rotq_hifi4 \
         imgfast_gs8_rotq16_hifi4 \
         imgfast_gs16_crop_hifi4 \
         imgfast_gs16_rot0_16_hifi4 \
         imgfast_gs16_rot0_hifi4 \
         imgfast_gs16_rot180_16_hifi4 \
         imgfast_gs16_rot180_hifi4 \
         imgfast_gs16_rot270_16_hifi4 \
         imgfast_gs16_rot270_hifi4 \
         imgfast_gs16_rot90_16_hifi4 \
         imgfast_gs16_rot90_hifi4 \
         imgfast_gs16_rotq_hifi4 \
         imgfast_gs16_rotq16_hifi4 \
         imgfastrotate_gu8_hifi4 \
         imgfastrotate_gs8_hifi4 \
         imgfastrotate_gs16_hifi4 \
         imgrotate_getCoef_hifi4 \
         imgrotate_gu8_hifi4 \
         imgrotate_gs8_hifi4 \
         imgrotate_gs16_hifi4 \
         imgrotate_process_hifi4 \
         cmatcholtransformR_hifi4 \
         cmatcholbkwsubstnxn_hifi4 \
         cmatcholfwdsubstnxn_hifi4 \
         cmatcholfwdsubstrec_hifi4 \
         cmatcholdiagUpd_hifi4 \
         cmatcholdecomp10x10_32x32_hifi4 \
         cmatcholdecomp4x4_32x32_hifi4 \
         cmatcholdecomp6x6_32x32_hifi4 \
         cmatcholdecomp8x8_32x32_hifi4 \
         cmatcholbkwsubst10x10_32x32_hifi4 \
         cmatcholbkwsubst4x4_32x32_hifi4 \
         cmatcholbkwsubst6x6_32x32_hifi4 \
         cmatcholbkwsubst8x8_32x32_hifi4 \
         cmatcholfwdsubst10x10_32x32_hifi4 \
         cmatcholfwdsubst4x4_32x32_hifi4 \
         cmatcholfwdsubst6x6_32x32_hifi4 \
         cmatcholfwdsubst8x8_32x32_hifi4 \
         cmatcholmmsesolver10x10_32x32_hifi4 \
         cmatcholmmsesolver4x4_32x32_hifi4 \
         cmatcholmmsesolver6x6_32x32_hifi4 \
         cmatcholmmsesolver8x8_32x32_hifi4 \
         cmatcholpseudoinv10x10_32x32_hifi4 \
         cmatcholpseudoinv4x4_32x32_hifi4 \
         cmatcholpseudoinv6x6_32x32_hifi4 \
         cmatcholpseudoinv8x8_32x32_hifi4 \
         cmatcholpreprocess10x10_32x32_hifi4 \
         cmatcholpreprocess4x4_32x32_hifi4 \
         cmatcholpreprocess6x6_32x32_hifi4 \
         cmatcholpreprocess8x8_32x32_hifi4 \
         matcholbkwsubstnxn_hifi4 \
         matcholfwdsubstnxn_hifi4 \
         matcholfwdsubstrec_hifi4 \
         matcholdiagUpd_hifi4 \
         matcholdecomp10x10_32x32_hifi4 \
         matcholdecomp4x4_32x32_hifi4 \
         matcholdecomp6x6_32x32_hifi4 \
         matcholdecomp8x8_32x32_hifi4 \
         matcholbkwsubst10x10_32x32_hifi4 \
         matcholbkwsubst4x4_32x32_hifi4 \
         matcholbkwsubst6x6_32x32_hifi4 \
         matcholbkwsubst8x8_32x32_hifi4 \
         matcholfwdsubst10x10_32x32_hifi4 \
         matcholfwdsubst4x4_32x32_hifi4 \
         matcholfwdsubst6x6_32x32_hifi4 \
         matcholfwdsubst8x8_32x32_hifi4 \
         matcholmmsesolver10x10_32x32_hifi4 \
         matcholmmsesolver4x4_32x32_hifi4 \
         matcholmmsesolver6x6_32x32_hifi4 \
         matcholmmsesolver8x8_32x32_hifi4 \
         matcholpseudoinv10x10_32x32_hifi4 \
         matcholpseudoinv4x4_32x32_hifi4 \
         matcholpseudoinv6x6_32x32_hifi4 \
         matcholpseudoinv8x8_32x32_hifi4 \
         matcholpreprocess10x10_32x32_hifi4 \
         matcholpreprocess4x4_32x32_hifi4 \
         matcholpreprocess6x6_32x32_hifi4 \
         matcholpreprocess8x8_32x32_hifi4 \
         cmatcholdecomp10x10f_hifi4 \
         cmatcholdecomp4x4f_hifi4 \
         cmatcholdecomp6x6f_hifi4 \
         cmatcholdecomp8x8f_hifi4 \
         cmatcholbkwsubst10x10f_hifi4 \
         cmatcholbkwsubst4x4f_hifi4 \
         cmatcholbkwsubst6x6f_hifi4 \
         cmatcholbkwsubst8x8f_hifi4 \
         cmatcholfwdsubst10x10f_hifi4 \
         cmatcholfwdsubst4x4f_hifi4 \
         cmatcholfwdsubst6x6f_hifi4 \
         cmatcholfwdsubst8x8f_hifi4 \
         cmatcholmmsesolver10x10f_hifi4 \
         cmatcholmmsesolver4x4f_hifi4 \
         cmatcholmmsesolver6x6f_hifi4 \
         cmatcholmmsesolver8x8f_hifi4 \
         cmatcholpseudoinv10x10f_hifi4 \
         cmatcholpseudoinv4x4f_hifi4 \
         cmatcholpseudoinv6x6f_hifi4 \
         cmatcholpseudoinv8x8f_hifi4 \
         cmatcholpreprocess10x10f_hifi4 \
         cmatcholpreprocess4x4f_hifi4 \
         cmatcholpreprocess6x6f_hifi4 \
         cmatcholpreprocess8x8f_hifi4 \
         cplxcholDiagUpdf_hifi4 \
         cplxcholFwdrec10f_hifi4 \
         cplxcholFwdrec1f_hifi4 \
         cplxcholFwdrec2f_hifi4 \
         cplxcholFwdrec3f_hifi4 \
         cplxcholFwdrec4f_hifi4 \
         cplxcholFwdrec5f_hifi4 \
         cplxcholFwdrec6f_hifi4 \
         cplxcholFwdrec7f_hifi4 \
         cplxcholFwdrec8f_hifi4 \
         cplxcholFwdrec9f_hifi4 \
         matcholdecomp10x10f_hifi4 \
         matcholdecomp4x4f_hifi4 \
         matcholdecomp6x6f_hifi4 \
         matcholdecomp8x8f_hifi4 \
         matcholbkwsubst10x10f_hifi4 \
         matcholbkwsubst4x4f_hifi4 \
         matcholbkwsubst6x6f_hifi4 \
         matcholbkwsubst8x8f_hifi4 \
         matcholfwdsubst10x10f_hifi4 \
         matcholfwdsubst4x4f_hifi4 \
         matcholfwdsubst6x6f_hifi4 \
         matcholfwdsubst8x8f_hifi4 \
         matcholmmsesolver10x10f_hifi4 \
         matcholmmsesolver4x4f_hifi4 \
         matcholmmsesolver6x6f_hifi4 \
         matcholmmsesolver8x8f_hifi4 \
         matcholpseudoinv10x10f_hifi4 \
         matcholpseudoinv4x4f_hifi4 \
         matcholpseudoinv6x6f_hifi4 \
         matcholpseudoinv8x8f_hifi4 \
         matcholpreprocess10x10f_hifi4 \
         matcholpreprocess4x4f_hifi4 \
         matcholpreprocess6x6f_hifi4 \
         matcholpreprocess8x8f_hifi4 \
         realcholDiagUpdf_hifi4 \
         realcholFwdrec10f_hifi4 \
         realcholFwdrec1f_hifi4 \
         realcholFwdrec2f_hifi4 \
         realcholFwdrec3f_hifi4 \
         realcholFwdrec4f_hifi4 \
         realcholFwdrec5f_hifi4 \
         realcholFwdrec6f_hifi4 \
         realcholFwdrec7f_hifi4 \
         realcholFwdrec8f_hifi4 \
         realcholFwdrec9f_hifi4 \
         NatureDSP_Signal_complex_id \
         NatureDSP_Signal_fft_id \
         NatureDSP_Signal_fir_id \
         NatureDSP_Signal_fit_id \
         NatureDSP_Signal_iir_id \
         NatureDSP_Signal_math_id \
         NatureDSP_Signal_matinv_id \
         NatureDSP_Signal_matop_id \
         NatureDSP_Signal_vector_id \
         NatureDSP_Signal_mfcc_id \
         NatureDSP_Signal_img_id \
         NatureDSP_Signal_fe \
         NatureDSP_Signal_isa_opt \
         feature \
         math_stdlib \
		logmel32x32_hifi4 \
		logmel32x32_process_hifi4 \
		logmel_common_hifi4 \
		logmelf_hifi4 \
		logmelf_process_hifi4 \
		mfcc32x32_compDctMatrix_hifi4 \
		mfcc32x32_compLifterCoefs_hifi4 \
		mfcc32x32_hifi4 \
		mfcc32x32_preemph_hifi4 \
		mfcc32x32_process_hifi4 \
		mfcc32x32_remdc_hifi4 \
		mfcc32x32_tbl \
		mfcc32x32_vecmpy_hifi4 \
		mfcc_common_hifi4 \
		mfccf_compDctMatrix_hifi4 \
		mfccf_compLifterCoefs_hifi4 \
		mfccf_hifi4 \
		mfccf_preemph_hifi4 \
		mfccf_process_hifi4 \
		mfccf_remdc_hifi4 \
		mfccf_vecmpy_hifi4 \
         version 

ifneq (0, $(HAVE_XTENSA))
#----------------------------------------------------------
# CFLAGS_FILE_module_name   += some compiler switches
# LLVMFLAGS_FILE_module_name+= some compiler switches
# select specific compiler switches for xcc and llvm(clang)
# compilers respectively
#----------------------------------------------------------
#CFLAGS_FILE_module_name   += some compiler switches
#CFLAGS_FILE_bkfir24x24p_process_hifi4 +=-O0


LLVMFLAGS_FILE_vec_shiftf_hifi4          += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_antilog10f_hifi4      += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_antilog2f_hifi4       += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_antilognf_hifi4       += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_atan2f_hifi4          += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_atanf_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_cosinef_hifi4         += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_float2int_hifi4       += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_int2float_hifi4       += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_log10f_hifi4          += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_log2f_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_lognf_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_sinef_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_scl_tanf_hifi4            += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_alog10f_hifi4         += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_alog2f_hifi4          += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_alognf_hifi4          += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_atan2f_hifi4          += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_atanf_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_cosinef_hifi4         += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_float2int_hifi4       += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_int2float_hifi4       += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_log10f_hifi4          += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_log2f_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_lognf_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_sinef_hifi4           += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_vec_tanf_hifi4            += -fno-unsafe-math-optimizations

LLVMFLAGS_FILE_latrf_hifi4               += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf1_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf2_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf3_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf4_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf5_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf6_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf7_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrf8_hifi4              += -fno-unsafe-math-optimizations
LLVMFLAGS_FILE_latrfX_hifi4              += -fno-unsafe-math-optimizations

endif


			
default:	all

.PHONY:	all clean cleanest cleanlib debug ferret default profile release
.PHONY: $(foreach proj,$(EXTERNAL_PROJECTS), $(proj))
.PHONY: $(foreach proj,$(EXTERNAL_PROJECTS), $(proj)/clean)
.PHONY: $(foreach proj,$(EXTERNAL_PROJECTS), $(proj)/cleanest)

debug:		all
ferret:		all
profile:	all
release:	all 

ifeq (0,$(CPP)) ######################################
ifeq (LLVM,$(LANG))
PLATFORM		:= $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)_llvm
else
PLATFORM		:= $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)
endif
else #################################################
PLATFORM		:= $(XTCORE)$(MEM_SUFFIX)-C++
endif ################################################

ifneq (, $(findstring cstub, $(PLATFORM))) ###########
SOLID_LIBRARY		:= 1
CCFLAGS_DIR_fft_files	+= -Wno-uninitialized -Wno-strict-aliasing
CCFLAGS_DIR_cstub 	+= -Wno-uninitialized -Wno-strict-aliasing
CCFLAGS_DIR_library	+= -Wno-strict-aliasing
CCFLAGS_DIR_P2_files	+= -Wno-strict-aliasing
else #################################################
SOLID_LIBRARY		:= 0
CCFLAGS_DIR_fft_files += -OPT:alias=restrict
CCFLAGS_DIR_P2_files += -LNO:simd
endif ################################################
###############################################################
##                      Host OS detection                    ##
###############################################################

WINDIR+=$(windir)

ifeq (,$(WINDIR))
# Unix-like system
ULIS=1
override HOSTOS=$(shell uname)
else
ULIS=0
override HOSTOS=win32
SHELL=cmd
endif

###############################################################
##              Platform-specific configuration              ##
###############################################################
CSTUB_MODULES += cstub-Xm_$(XTCORE)

ifeq (,$(XTENSA_SYSTEM))
CSTUB_INCDIRS	+= $(CSTUB_DIR) $(CSTUB_DIR)/xtensa/tie $(CSTUB_DIR)/xtensa/config
CSTUB_SRCDIRS	+= $(CSTUB_DIR)
CFLAGS_FILE_$(subst -,_,$(CSTUB_MODULES)) = -Dasm=__asm__ -O0
LLVMFLAGS_FILE_$(subst -,_,$(CSTUB_MODULES)) = -Dasm=__asm__ -O0
else
CSTUB_INCDIRS += $(XTENSA_SYSTEM)/../src/cstub $(XTENSA_SYSTEM)/../src/cstub/xtensa/tie $(XTENSA_SYSTEM)/../xtensa-elf/arch/include/xtensa/config
CSTUB_SRCDIRS += $(XTENSA_SYSTEM)/../src/cstub
endif


ifneq (0, $(HAVE_XTENSA)) ############################
#############################################################
VARPREFIX.$(PLATFORM)		     = $(PLATFORM)
ifeq (LLVM,$(LANG))
CC.$(PLATFORM)      		     = xt-clang
LD.$(PLATFORM)      		     = xt-clang
else
CC.$(PLATFORM)      		     = xt-xcc 
LD.$(PLATFORM)      		     = xt-xcc
endif
AR.$(PLATFORM)      		     = xt-ar
#MAKE.$(PLATFORM)      		   = xt-make
OBJCOPY.$(PLATFORM)          = xt-objcopy
RUN.$(PLATFORM)     		     = xt-run
PROF.$(PLATFORM)      		   = xt-prof
CFLAGS.$(PLATFORM) 		      += -pipe 
CFLAGS.$(PLATFORM)         	+= 
CFLAGS.profile.$(PLATFORM)	+= -O3 
CFLAGS.release.$(PLATFORM)	+= -O3  
CFLAGS.ferret.$(PLATFORM)	+= -O3 -g 
ifeq (LLVM,$(LANG))
CFLAGS.profile.$(PLATFORM)	+= -mlongcalls -fno-builtin
CFLAGS.release.$(PLATFORM)	+= -mlongcalls -fno-builtin
CFLAGS.ferret.$(PLATFORM)	+= -mlongcalls -fno-builtin
else
CFLAGS.profile.$(PLATFORM)	+= -LNO:simd -mlongcalls -fno-builtin
CFLAGS.release.$(PLATFORM)	+= -LNO:simd -mlongcalls -fno-builtin
CFLAGS.ferret.$(PLATFORM)	+= -LNO:simd -mlongcalls -fno-builtin
endif

CFLAGS.$(PLATFORM)	+= $(CFLAGS_EXT)
RUNFLAGS.debug.$(PLATFORM)	 = --turbo
RUNFLAGS.ferret.$(PLATFORM)	 = --turbo
RUNFLAGS.release.$(PLATFORM) = --turbo
TARGETOS.$(PLATFORM)		     = Xtensa
else
VARPREFIX.$(PLATFORM)		 = $(PLATFORM)
ifeq (1, $(CPP))
CC.$(PLATFORM)      		 = g++
LD.$(PLATFORM)      		 = g++
else
VARPREFIX.$(PLATFORM)		 = $(PLATFORM)
CC.$(PLATFORM)      		 = gcc
LD.$(PLATFORM)      		 = gcc
endif
AR.$(PLATFORM)      		 = ar
#MAKE.$(PLATFORM)      		 = make
OBJCOPY.$(PLATFORM)              = objcopy
RUN.$(PLATFORM)     		 = 
PROF.$(PLATFORM)      		 = 
CFLAGS.$(PLATFORM)         	 += -m32 -pipe 
CFLAGS.$(PLATFORM) 		 += -DXTENSA_CORE=Xm_$(XTCORE)
CFLAGS.$(PLATFORM)		+= -D__USE_XOPEN2K8	# patch for bug in gcc
ifeq (XRC_B64_LP,$(XTCORE))
CFLAGS.$(PLATFORM)         	 += -DB64_LP
endif
INCDIRS			 	 += $(CSTUB_INCDIRS)
#SRCDIRS			 += $(CSTUB_SRCDIRS)
ifeq (1, $(COMPILE_CSTUB))
MODULES				 += $(CSTUB_MODULES)
endif
CFLAGS.profile.$(PLATFORM)       += -O2 
CFLAGS.release.$(PLATFORM)	 += -O2 
RUNFLAGS.debug.$(PLATFORM)	 = 
RUNFLAGS.release.$(PLATFORM)	 = 
RUNFLAGS.ferret.$(PLATFORM)	 = 
TARGETOS.$(PLATFORM)		 = Linux
LFLAGS.$(PLATFORM) 		 += -m32
endif

###############################################################
##                     Well-known suffixes                   ##
###############################################################
#POSTFIX.$(PLATFORM)-Xtensa-profile.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-release.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-debug.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-ferret.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-profile.1=_$(PLATFORM)

###############################################################
##                 Host-specific configuration               ##
###############################################################

ifeq (1,$(ULIS)) ##############################################
ECHO=echo
override RM=rm -f
CP=cp
PSEP=/
NULLDEV=/dev/null
MKDIR=mkdir
JOIN=&&
else ##########################################################
ECHO=echo
override RM=del
CP=copy
PSEP=\\#
NULLDEV=NUL
MKDIR=mkdir
JOIN=&
endif #########################################################

###############################################################
##                     Target OS definition                  ##
###############################################################

TARGETOS=$(firstword $(TARGETOS.$(PLATFORM)) $(HOSTOS))

VARPREFIX.Darwin=MAC
VARPREFIX.iOS=MAC
VARPREFIX.linux=LINUX
VARPREFIX.Linux=LINUX
VARPREFIX.unix-like.1=UNIX_LIKE
VARPREFIX.win32=WIN32

###############################################################
##                  Select build configuration               ##
###############################################################

CONFIG=$(firstword $(findstring release, $(MAKECMDGOALS)) $(findstring debug,\
 $(MAKECMDGOALS)) $(findstring ferret, $(MAKECMDGOALS)) $(findstring profile, $(MAKECMDGOALS)) release)

CONFIG_ID	= $(PLATFORM)-$(TARGETOS)-$(CONFIG)
DLIBS_CONFIG_ID	= $(PLATFORM)-$(TARGETOS)-$(firstword \
		  $(DLIBS_CONFIG) $(CONFIG))

# Intermediate directory

INTERMEDIATE_DIR= $(CONFIG_ID)
USE_SUFFIXES	= $(firstword $(INTEGRIT_SUFFIXES) $(INTEGRIT_OLD_SUFFIXES))
POSTFIX		= $(firstword \
	$(POSTFIX.$(CONFIG_ID).$(USE_SUFFIXES)) -$(CONFIG_ID))
POSTFIX_DLC	= $(firstword \
	$(POSTFIX.$(DLIBS_CONFIG_ID).$(USE_SUFFIXES)) -$(DLIBS_CONFIG_ID))

$(INTERMEDIATE_DIR)/dirinfo.txt:
	-$(MKDIR) $(INTERMEDIATE_DIR)
	-$(ECHO) Intermediate directory.>$(INTERMEDIATE_DIR)/dirinfo.txt

###############################################################
##                     Configure DevTools                    ##
###############################################################
CC			= $(firstword $(CC.$(PLATFORM)) gcc)
LD			= $(firstword $(LD.$(PLATFORM)) gcc)
AR			= $(firstword $(AR.$(PLATFORM)) ar)
MAKE			= $(firstword $(MAKE.$(PLATFORM)) make)
OBJCOPY			= $(firstword $(OBJCOPY.$(PLATFORM)) objcopy)
RUN			= $(RUN.$(PLATFORM))
PROF			= $(firstword $(PROF.$(PLATFORM)) gprof)

CFLAGS_DEPENDENT	+= $(CFLAGS.$(PLATFORM))

CFLAGS.default.debug.	= $(CFLAGS.default.debug) -g -O0 -mtext-section-literals -mlongcalls -D_DEBUG 
CFLAGS.default.ferret.	= $(CFLAGS.default.ferret) -g -O0 -mtext-section-literals -mlongcalls -D_DEBUG 
CFLAGS.default.profile.	= $(CFLAGS.default.profile) -g -pg -O3
CFLAGS.default.release.	= $(CFLAGS.default.release) -O3
LFLAGS.default.debug.	= $(LFLAGS.default.debug)
LFLAGS.default.ferret.	= $(LFLAGS.default.ferret)
LFLAGS.default.profile.	= $(LFLAGS.default.profile) -g -pg -O
LFLAGS.default.release.EXECUTABLE = -s
LFLAGS.default.ferret.EXECUTABLE = 
LFLAGS.default.release.	= $(LFLAGS.default.release) -O \
			$(LFLAGS.default.release.$(BUILD_TYPE))
LFLAGS.default.ferret.	= $(LFLAGS.default.ferret) -O \
			$(LFLAGS.default.ferret.$(BUILD_TYPE))

CFLAGS.debug		= $(CFLAGS.debug.$(PLATFORM)) \
			$(CFLAGS.default.debug.$(CFLAGS.debug.$(PLATFORM)))
CFLAGS.ferret		= $(CFLAGS.ferret.$(PLATFORM)) \
			$(CFLAGS.default.ferret.$(CFLAGS.ferret.$(PLATFORM)))
CFLAGS.profile		= $(CFLAGS.profile.$(PLATFORM)) \
			$(CFLAGS.default.profile.$(CFLAGS.profile.$(PLATFORM)))
CFLAGS.release		= $(CFLAGS.release.$(PLATFORM)) \
			$(CFLAGS.default.release.$(CFLAGS.release.$(PLATFORM)))

CFLAGADDINCDIR		= $(firstword $(CFLAGADDINCDIR.$(PLATFORM)) -I)
CFLAGMAKEDEP		= $(firstword $(CFLAGMAKEDEP.$(PLATFORM)) -MM)
CFLAGMAKEOBJ		= $(firstword $(CFLAGMAKEDEP.$(PLATFORM)) -c)
#CFLAGSAVETEMP.1		= $(firstword $(CFLAGSAVETEMP.$(PLATFORM)) -save-temps)
CFLAGSAVETEMP.1		=  -save-temps -g

LFLAGS_DEPENDENT	+= $(LFLAGS.$(PLATFORM))
LFLAGS.debug		= $(LFLAGS.debug.$(PLATFORM)) \
			$(LFLAGS.default.debug.$(LFLAGS.debug.$(PLATFORM)))
LFLAGS.debug		= $(LFLAGS.ferret.$(PLATFORM)) \
			$(LFLAGS.default.ferret.$(LFLAGS.ferret.$(PLATFORM)))
LFLAGS.profile		= $(LFLAGS.profile.$(PLATFORM)) \
			$(LFLAGS.default.profile.$(LFLAGS.profile.$(PLATFORM)))
LFLAGS.release		= $(LFLAGS.release.$(PLATFORM)) \
			$(LFLAGS.default.release.$(LFLAGS.release.$(PLATFORM)))

LFLAGS.LIBRARY.1	+=--no-standard-libraries
LFLAGS.LIBRARY.1	+=-Wl,-r,--retain-symbols-file,$(TARGET_SYMFILE)

ARFLAGS			= $(firstword $(ARFLAGS.$(PLATFORM)) rcs)

RUNFLAGS		+= $(RUNFLAGS.$(PLATFORM))
RUNFLAGS		+= $(RUNFLAGS.$(CONFIG).$(PLATFORM))

PROFFLAGS		+= $(PROFFLAGS.$(PLATFORM))

LFLAGADDLIBDIR		= $(firstword $(LFLAGADDLIBDIR.$(PLATFORM)) -L)

INCDIRS			+= $($(VARPREFIX.$(PLATFORM))_INCDIRS)
SRCDIRS			+= $($(VARPREFIX.$(PLATFORM))_SRCDIRS)
LIBDIRS			+= $($(VARPREFIX.$(PLATFORM))_LIBDIRS)
SRCDIRS			+= $(SRCDIR)

# Header search options
CFLAGS_DEPENDENT += $(foreach idir, $(INCDIRS), $(CFLAGADDINCDIR) $(idir))

CFLAGS_DEPENDENT	+= $(CFLAGSAVETEMP.$(SAVETEMPS))

# Library search options
LFLAGS_DEPENDENT += $(foreach ldir, $(LIBDIRS), $(LFLAGADDLIBDIR) $(ldir))

# Configuration-dependent options
RFLAGS.release = -D__RENAMING__
#RFLAGS.debug = -D__RENAMING__
RFLAGS.ferret =
RFLAGS = $(RFLAGS.$(CONFIG))

CFLAGS_DEPENDENT += $(CFLAGS.$(CONFIG))
LFLAGS_DEPENDENT += $(LFLAGS.$(CONFIG)) $(LFLAGS.$(BUILD_TYPE).$(SOLID_LIBRARY))
CFLAGS_DEPENDENT += -Wall -Wno-return-type
ifeq (0, $(CPP))
CFLAGS_DEPENDENT += -std=c99
CFLAGS_DEPENDENT += -Werror
else
CFLAGS_DEPENDENT += -Wno-error
#CFLAGS_DEPENDENT += -Werror
endif
CFLAGS_DEPENDENT += $(RFLAGS)

# support memory modeling and data cache
ifeq (2,$(MEM_MODEL))
CFLAGS.$(PLATFORM)		+= -DDATA_CACHE=$(DATA_CACHE) -DMEM_MODEL
else
ifeq (1,$(MEM_MODEL)) 
CFLAGS.$(PLATFORM)		+= -DMEM_MODEL -DDATA_CACHE=0
else
CFLAGS.$(PLATFORM)		+= -DDATA_CACHE=0
endif
endif
###############################################################
##                     Create target spec                    ##
###############################################################

TARGET_SUFFIX.win32	= .exe

TARGET_SUFFIX		= $(firstword $(TARGET_SUFFIX.$(BUILD_TYPE))\
			  $(TARGET_SUFFIX.$(HOSTOS)))

TARGET_SPEC.LIBRARY	= $(OUTDIR)/$(TARGET)$(POSTFIX).a
TARGET_SPEC.EXECUTABLE	= $(OUTDIR)/$(TARGET)$(POSTFIX)$(TARGET_SUFFIX)

TARGET_INTERMEDIATE.LIBRARY.1.1	= $(INTERMEDIATE_DIR)/$(TARGET).solid.1.o
TARGET_INTERMEDIATE.LIBRARY.1.2	= $(INTERMEDIATE_DIR)/$(TARGET).solid.2.o
TARGET_INTERMEDIATE.LIBRARY.1.3	= $(INTERMEDIATE_DIR)/$(TARGET).solid.o

###############################################################
##                    Complete module list                   ##
###############################################################

MODULES += $($(VARPREFIX.$(PLATFORM))_MODULES)
MODULES += $($(VARPREFIX.$(TARGETOS))_MODULES)
MODULES += $($(VARPREFIX.unix-like.$(ULIS))_MODULES)
# Keep compatibility with IntegrIT makefiles below 0.30
MODULES += $(OBJS:.o=)
MODULES += $($(VARPREFIX.$(PLATFORM))_OBJS:.o=)
MODULES += $($(VARPREFIX.$(TARGETOS))_OBJS:.o=)
MODULES += $($(VARPREFIX.unix-like.$(ULIS))_OBJS:.o=)

#     WIN32_MODULES, UNIX_LIKE_MODULES, MAC_MODULES

###############################################################
##           Generate list of objects to be built            ##
###############################################################

ALL_OBJECTS = $(foreach module, $(MODULES), $(INTERMEDIATE_DIR)/$(module).o)

###############################################################
##        Generate list of extra libraries to be linked      ##
###############################################################

ALL_DLIBS=$(foreach lib,$(DLIBS),$(lib)$(POSTFIX_DLC).a)

###############################################################
##                     Dependency generator                  ##
###############################################################

ifneq (1,$(CPP))
vpath %.c $(SRCDIRS) $(CSTUB_SRCDIRS)
else
vpath %.c $(SRCDIRS)
vpath %.cpp $(SRCDIRS) $(CSTUB_SRCDIRS)
endif

COMMAND_DEPGEN =  $(CC)
COMMAND_DEPGEN += $(CFLAGS_DEPENDENT)
COMMAND_DEPGEN += $(CFLAGS_DIR_$(subst -,_,$(notdir $(dir $(1)))))
ifeq (LLVM,$(LANG))
COMMAND_DEPGEN += $(LLVMFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
else
COMMAND_DEPGEN += $(CFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
endif
COMMAND_DEPGEN += $(CFLAGS)
COMMAND_DEPGEN += $(CFLAGMAKEDEP)
COMMAND_DEPGEN += -MQ $(2)
COMMAND_DEPGEN += -MQ $(2:.d=.o)
COMMAND_DEPGEN += $(1)
COMMAND_DEPGEN += -o $(2)

$(INTERMEDIATE_DIR)/%.d: %.c $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_DEPGEN,$<,$@)

$(INTERMEDIATE_DIR)/%.d: %.cpp $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_DEPGEN,$<,$@)

# Include dependencies
-include $(addprefix $(INTERMEDIATE_DIR)/, $(addsuffix .d, $(MODULES)))

###############################################################
##                        Compiling rules                    ##
###############################################################

COMMAND_COMPILE =  $(CC)
COMMAND_COMPILE += $(CFLAGS_DEPENDENT)
COMMAND_COMPILE += $(CFLAGS_DIR_$(subst -,_,$(notdir $(dir $(1)))))
ifeq (LLVM,$(LANG))
COMMAND_COMPILE += $(LLVMFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
else
COMMAND_COMPILE += $(CFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
endif
COMMAND_COMPILE += $(CFLAGS)
COMMAND_COMPILE += $(CFLAGMAKEOBJ)
COMMAND_COMPILE += $(1)
COMMAND_COMPILE += -o $(2) 

$(INTERMEDIATE_DIR)/%.o: %.c $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_COMPILE,$<,$@)

$(INTERMEDIATE_DIR)/%.o: %.cpp $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_COMPILE,$<,$@)

###############################################################
##                    Linking/archiving rules                ##
###############################################################

ALL_LIBS = $(ALL_DLIBS) $(foreach lib,$(LIBS),-l$(lib))

LD_LIST.  = $(ALL_OBJECTS)
LD_LIST.0 = $(ALL_OBJECTS)
LD_LIST.1 = $(TARGET_INTERMEDIATE.$(BUILD_TYPE).$(SOLID_LIBRARY).3)

$(TARGET_INTERMEDIATE.LIBRARY.1.1): $(ALL_OBJECTS) $(ALL_DLIBS)
	$(LD) $(LFLAGS_DEPENDENT) $(LFLAGS) -o \
	$(TARGET_INTERMEDIATE.LIBRARY.1.1) $(ALL_OBJECTS) $(ALL_LIBS)

$(TARGET_INTERMEDIATE.LIBRARY.1.2):	$(TARGET_INTERMEDIATE.LIBRARY.1.1)
	$(OBJCOPY) --keep-global-symbols=$(TARGET_SYMFILE) \
	$(TARGET_INTERMEDIATE.LIBRARY.1.1) $(TARGET_INTERMEDIATE.LIBRARY.1.2)

$(TARGET_INTERMEDIATE.LIBRARY.1.3):	$(TARGET_INTERMEDIATE.LIBRARY.1.2)
	$(OBJCOPY) --strip-unneeded \
	$(TARGET_INTERMEDIATE.LIBRARY.1.2) $(TARGET_INTERMEDIATE.LIBRARY.1.3)

genlist.cmd:
	$(ECHO) @ECHO OFF>genlist.cmd
	$(ECHO) IF EXIST %%2 DEL %%2>>genlist.cmd
	$(ECHO) FOR %%%%i IN (%%1\*.o) DO (>>genlist.cmd
	$(ECHO)   ECHO %%%%i^>^>%%2>>genlist.cmd
	$(ECHO) )>>genlist.cmd

GENLIST = $(if $(WINDIR),genlist.cmd,)

$(INTERMEDIATE_DIR)/$(TARGET).ld: $(LD_LIST.$(SOLID_LIBRARY)) $(GENLIST)
ifeq (,$(WINDIR))
	$(shell $(RM) $@) $(foreach mod,$(LD_LIST.$(SOLID_LIBRARY)),$(shell $(ECHO) $(mod)>>$@))
else
	$(GENLIST) $(INTERMEDIATE_DIR) $(subst /,$(PSEP),$@)
endif

$(TARGET_SPEC.LIBRARY):	$(EXTERNAL_PROJECTS) $(LD_LIST.$(SOLID_LIBRARY)) $(INTERMEDIATE_DIR)/$(TARGET).ld
	$(RM) $(subst /,$(PSEP),$@)
	$(AR) rcs $@ @$(INTERMEDIATE_DIR)/$(TARGET).ld

$(TARGET_SPEC.EXECUTABLE): $(EXTERNAL_PROJECTS) $(ALL_OBJECTS) $(ALL_DLIBS)
	$(LD) $(LFLAGS_DEPENDENT) $(LFLAGS) -o \
	$(TARGET_SPEC.EXECUTABLE) $(ALL_OBJECTS) $(ALL_LIBS)

all:	$(TARGET_SPEC.$(BUILD_TYPE))

$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)):
	make -C $@ CONFIG=$(CONFIG) all

###############################################################
##                         Clean rules                       ##
###############################################################
$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/clean):
	make -C $(dir $@) CONFIG=$(CONFIG) clean
$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/cleanest):
	make -C $(dir $@) CONFIG=$(CONFIG) cleanest

clean:	$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/clean)
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.o)
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.$(BUILD_TYPE)))

cleanest:	$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/cleanest)
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.d)
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.o)
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.$(BUILD_TYPE)))

cleanlib:
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.LIBRARY)) 2>$(NULLDEV)
