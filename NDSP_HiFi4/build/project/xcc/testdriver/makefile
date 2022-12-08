# makefile for building IntegrIT(R) DSP Library
# Version 0.34 tuned for Cadence HiFi cores
# ---------------------------------------------------
# Target platform: gcc/xcc
#
#  make invocation:
#  make clean   - clean up all temporary files
#  make depends - refresh dependencies (makedeps.mk)
#  make all     - make depends and release
#-----------------------------------------------------
#  build under xcc
#  make sure that all neccesary Xtensa system variables 
#  and pathes configured properly. makefile uses 
#  XTENSA_CORE system variable to detect the name of core
#  so, if makefile is running from the command-line shell
#  executed from proper processor config of Xtensa Exmplorer 
#  it does not require setting of any additional parameters
#
#  System variables XTENSA_CORE,XTENSA_SYSTEM should.
#  be properly configured for correct procedure.
#  if XTENSA_SYSTEM is undefined, g++ build procedure
#  will use cstub files and headers from local
#  directory, not from system path
#
#  Other variables:
#  LANG=C++    build under g++, otherwise xcc used
#  LANG=LLVM   build under llvm Xtensa compiler
#  SAVETEMPS=1 save temporary files
#  MEM_MODEL=1,2 enables memory modeling
#-----------------------------------------------------  
#  build under gcc require unsetting of variable XTENSA_SYSTEM
#  (either by command line unset XTENSA_SYSTEM or via make 
#  command line option -e XTENSA_SYSTEM=)
#  and selecting proper core by XTENSA_CORE variable. 
#  By default, it use C language, C++ will be selected by 
#  using LANG=C++ system variable.
# 
#  examples:
#  1. build with with cstubs :
#  make -e XTENSA_SYSTEM= XTENSA_CORE=XRC_B64EP_PM_20apr
#  2. build under g++ with cstubs:
#  make -e XTENSA_SYSTEM= XTENSA_CORE=XRC_B64EP_PM_20apr LANG=C++
#-----------------------------------------------------

SRCDIRS  =	../../../../testdriver \
             ../../../../testdriver/env \
             ../../../../testdriver/pkg/audio \
             ../../../../testdriver/pkg/complex \
             ../../../../testdriver/pkg/dct \
             ../../../../testdriver/pkg/fft \
             ../../../../testdriver/pkg/fir \
             ../../../../testdriver/pkg/chol \
             ../../../../testdriver/pkg/fit \
             ../../../../testdriver/pkg/iir \
             ../../../../testdriver/pkg/img \
             ../../../../testdriver/pkg/math \
             ../../../../testdriver/pkg/matinv \
             ../../../../testdriver/pkg/matop \
             ../../../../testdriver/pkg/vector

OUTDIR  =	../../../bin

INCDIRS =	../../../../library/include		\
		../../../../library/include_private	\
		      ../../../../testdriver/include

TARGET := testdriver
LIBS   := m
DLIBS  := ./../../../../build/bin/NatureDSP_Signal 

BUILD_TYPE		:= EXECUTABLE
INTEGRIT_OLD_SUFFIXES	:= 1
COMPILE_CSTUB		:= 0
LSP			:=
EXTERNAL_PROJECTS = ./../library

# check if c++ compilation required
ifeq (C++,$(LANG))
CPP = 1
HAVE_XTENSA = 0
else
CPP = 0
HAVE_XTENSA = 1
endif  

ifeq (,$(XTENSA_SYSTEM))
XTCORE := $(notdir $(wildcard ../../../../cstub/*-0.cpp))
XTCORE := $(subst -0.cpp,,$(XTCORE))
XTCORE := $(subst cstub-Xm_,,$(XTCORE))
else
XTCORE := $(XTENSA_CORE)
endif

EXTERNAL_PROJECTS.0 = ../cstub
EXTERNAL_PROJECTS += $(EXTERNAL_PROJECTS.$(HAVE_XTENSA))

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


MODULES = \
          mips \
          main \
          addr2name \
          float16 \
          fpstat \
          malloc16 \
          profiler \
          rms \
          testcase \
          testeng \
          testeng_errh \
          testeng_load_fxn \
          testeng_process_fxn \
          utils \
          vectools \
          vreport \
          mips_complexs \
          mips_complexv \
          test_complexs \
          test_complexv \
          mips_dct \
          test_dct \
          mips_chol \
          test_chol \
          mips_img \
          test_img \
          test_imgfft \
          testeng_dct \
          mips_conv2d \
          test_conv2d \
          mips_cfft \
          mips_cfftie \
          mips_cnfft \
          mips_rfft \
          mips_rfftie \
          mips_rnfft \
          test_cfft \
          test_rfft \
          testeng_fft \
          mips_firblk \
          mips_firdec \
          mips_firint \
          mips_firother \
          test_firblk \
          test_firdec \
          test_firint \
          test_firother \
          testeng_fir \
          testeng_fir_old \
          mips_pfit \
          test_pfit \
          mips_iirbq \
          mips_iirlt \
          test_iirbq \
          test_iirlt \
          testeng_iir_lat \
          testeng_iir_old \
          mips_maths \
          accuracy_maths \
          mips_mathv \
          test_maths \
          test_mathv \
          testeng_math \
          mips_matinv \
          test_matinv \
          testeng_matinv \
          mips_matop \
          test_matop \
          testeng_matop \
          mips_vector \
          test_vector \
          testeng_vector \
	  test_fft_spectrum \
	  mips_fft_spectrum \
	mips_mfcc \
	test_mfcc \
	testeng_mfcc \
	testeng_logmel \
	mips_ef \
	test_ef 

	
 
ifneq (0, $(HAVE_XTENSA))
#----------------------------------------------------------
# CFLAGS_FILE_module_name   += some compiler switches
# LLVMFLAGS_FILE_module_name+= some compiler switches
# select specific compiler switches for xcc and llvm(clang)
# compilers respectively
#----------------------------------------------------------
#CFLAGS_FILE_module_name   += some compiler switches

CFLAGS_FILE_mips_img += -mtext-section-literals
LLVMFLAGS_FILE_mips_img += -mtext-section-literals 

ifneq (,$(DRAM0))
DRAM0_MEM="$(DRAM0).data"
CFLAGS_FILE_mips += '-DDRAM0=$(DRAM0_MEM)'
endif
ifneq (,$(DRAM1))
DRAM1_MEM="$(DRAM1).data"
CFLAGS_FILE_mips += '-DDRAM1=$(DRAM1_MEM)'
endif
else
CFLAGS_FILE_testeng_errh += -DXTENSA_CORE=Xm_$(XTCORE)
endif

######### special definition for proper build under xws ##############
ifeq (0,0)
CFLAGS_FILE_mips += 
LLVMFLAGS_FILE_mips += 
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
PLATFORM     := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)_llvm
PLATFORM_LIB := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)_llvm
else
PLATFORM     := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)
PLATFORM_LIB := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)
endif
else #################################################
PLATFORM     := $(XTCORE)$(MEM_SUFFIX)-C++
PLATFORM_LIB := $(XTCORE)$(MEM_SUFFIX)-C++
endif ################################################

ifneq (, $(findstring cstub, $(PLATFORM))) ###########
SOLID_LIBRARY		:= 1
CCFLAGS_DIR_testdriver	+= -Wno-uninitialized
CCFLAGS_DIR_testdriver	+= -Wno-unused-result
CCFLAGS_DIR_cstub 		+= -Wno-uninitialized
CCFLAGS					+= -Wno-strict-aliasing
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
CSTUB_INCDIRS	+= ../../../../cstub
CSTUB_SRCDIRS	+= ../../../../cstub
ifeq (LLVM,$(LANG))
LLVMFLAGS_FILE_$(subst -,_,$(CSTUB_MODULES)) = -Dasm=__asm__ -O0
else
CFLAGS_FILE_$(subst -,_,$(CSTUB_MODULES)) = -Dasm=__asm__ -O0
endif
else
CSTUB_INCDIRS += $(XTENSA_SYSTEM)/../src/cstub $(XTENSA_SYSTEM)/../src/cstub/xtensa/tie
CSTUB_SRCDIRS += $(XTENSA_SYSTEM)/../src/cstub
endif

ifneq (0, $(HAVE_XTENSA)) ############################
##############################################################
VARPREFIX.$(PLATFORM)          = $(PLATFORM)
ifeq (LLVM,$(LANG))
CC.$(PLATFORM)      		     = xt-clang
LD.$(PLATFORM)      		     = xt-clang
else
CC.$(PLATFORM)      		     = xt-xcc 
LD.$(PLATFORM)      		     = xt-xcc
endif
AR.$(PLATFORM)		       = xt-ar
MAKE.$(PLATFORM)	       = xt-make
OBJCOPY.$(PLATFORM)	       = xt-objcopy
RUN.$(PLATFORM)		       = xt-run
RUNMEM.$(PLATFORM)	       = xt-run --mem_model
RUNFAST.$(PLATFORM)	       = xt-run --turbo
PROF.$(PLATFORM)	       = xt-prof
ARCH.$(PLATFORM)	       = $(PLATFORM)
CFLAGS.$(PLATFORM)	       += -fmessage-length=0 -mcoproc -mlongcalls  -fno-builtin
CFLAGS.profile.$(PLATFORM)     += -O3 
CFLAGS.release.$(PLATFORM)     += -O3 
ifeq (LLVM,$(LANG))
else
CFLAGS.$(PLATFORM)	       += -INLINE:requested
CFLAGS.profile.$(PLATFORM)     += -LNO:simd
CFLAGS.release.$(PLATFORM)     += -LNO:simd
endif
CFLAGS.$(PLATFORM)     += $(CFLAGS_EXT)
LFLAGS.profile.$(PLATFORM)     += --profile
RUNFLAGS.debug.$(PLATFORM)     = --turbo
RUNFLAGS.release.$(PLATFORM)   = --turbo
RUNFLAGS.ferret.$(PLATFORM)   = --turbo
TARGETOS.$(PLATFORM)	       = Xtensa
else
VARPREFIX.$(PLATFORM)		 = $(PLATFORM)
ifeq (1, $(CPP))
CC.$(PLATFORM)      		 = g++
LD.$(PLATFORM)      		 = g++
LFLAGS.$(PLATFORM) 		 += -Wl,--section-start=unused_section=0
CFLAGS.$(PLATFORM)		+=-D__USE_XOPEN2K8
else
VARPREFIX.$(PLATFORM)		 = $(PLATFORM)
CC.$(PLATFORM)      		 = gcc
LD.$(PLATFORM)      		 = gcc
LFLAGS.$(PLATFORM) 		 +=  -Wl,--section-start=unused_section=0
endif
AR.$(PLATFORM)      		 = ar
#MAKE.$(PLATFORM)      		 = make
OBJCOPY.$(PLATFORM)              = objcopy
RUN.$(PLATFORM)     		 = 
PROF.$(PLATFORM)      		 = 
CFLAGS.$(PLATFORM)         	 += -m32
CFLAGS.$(PLATFORM) 		 += -DXTENSA_CORE=Xm_$(XTCORE)
ifeq (XRC_B64_LP,$(XTCORE))
CFLAGS.$(PLATFORM)         	 += -DB64_LP
endif
INCDIRS			 	 += $(CSTUB_INCDIRS)
#SRCDIRS			 	 += $(CSTUB_SRCDIRS)
ifeq (1, $(COMPILE_CSTUB))
MODULES				 += $(CSTUB_MODULES)
endif
CFLAGS.profile.$(PLATFORM)       += -O3 
CFLAGS.release.$(PLATFORM)	 += -O3 
RUNFLAGS.debug.$(PLATFORM)	 = 
RUNFLAGS.release.$(PLATFORM)	 = 
RUNFLAGS.ferret.$(PLATFORM)	 = 
TARGETOS.$(PLATFORM)		 = Linux
LFLAGS.$(PLATFORM) 		 += -m32
DLIBS  += ../../../bin/cstub
endif

CFLAGS.release.$(PLATFORM)	 += -DUSE_ISA_PATCHES=$(USE_PATCHES)
###############################################################
##                     Well-known suffixes                   ##
###############################################################
#POSTFIX.$(PLATFORM)-Xtensa-release.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-debug.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-ferret.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-profile.1=_$(PLATFORM)

#POSTFIX_DLC.$(PLATFORM)-Xtensa-release.1=_$(PLATFORM_LIB)
#POSTFIX_DLC.$(PLATFORM)-Xtensa-debug.1=_$(PLATFORM_LIB)
#POSTFIX_DLC.$(PLATFORM)-Xtensa-profile.1=_$(PLATFORM_LIB)

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
ECHO=@echo
override RM=del
CP=copy
PSEP=\\\\#
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
DLIBS_CONFIG_ID	= $(PLATFORM_LIB)-$(TARGETOS)-$(firstword \
		  $(DLIBS_CONFIG) $(CONFIG))

# Intermediate directory

INTERMEDIATE_DIR= $(CONFIG_ID)
USE_SUFFIXES	= $(firstword $(INTEGRIT_SUFFIXES) $(INTEGRIT_OLD_SUFFIXES))
POSTFIX		= $(firstword \
	$(POSTFIX.$(CONFIG_ID).$(USE_SUFFIXES)) -$(CONFIG_ID))
POSTFIX_DLC	= $(firstword \
	$(POSTFIX_DLC.$(DLIBS_CONFIG_ID).$(USE_SUFFIXES)) \
	$(POSTFIX.$(DLIBS_CONFIG_ID).$(USE_SUFFIXES)) -$(DLIBS_CONFIG_ID))

$(INTERMEDIATE_DIR)/dirinfo.txt:
	-$(MKDIR) $(INTERMEDIATE_DIR)
	-$(ECHO) Intermediate directory.>$(INTERMEDIATE_DIR)/dirinfo.txt
	
###############################################################
##                Creating devtools options                  ##
###############################################################

IDIRS_OPT = $(foreach f,$(INCDIRS),-I $(f))
LDIRS_OPT = $(foreach f,$(LIBDIRS),-L $(f))

CC			= $(firstword $(CC.$(PLATFORM)) gcc)
LD			= $(firstword $(LD.$(PLATFORM)) gcc)
AR			= $(firstword $(AR.$(PLATFORM)) ar)
#MAKE			= $(firstword $(MAKE.$(PLATFORM)) make)
OBJCOPY			= $(firstword $(OBJCOPY.$(PLATFORM)) objcopy)
RUN			= $(RUN.$(PLATFORM))
PROF			= $(firstword $(PROF.$(PLATFORM)) gprof)

CFLAGS_DEPENDENT	+= $(CFLAGS.$(PLATFORM))

CFLAGS.default.debug.	= $(CFLAGS.default.debug) -g -O0 -D_DEBUG 
CFLAGS.default.profile.	= $(CFLAGS.default.profile) -g -pg -O3
CFLAGS.default.release.	= $(CFLAGS.default.release) -O3
CFLAGS.default.ferret.	= $(CFLAGS.default.ferret) -g -O0 -D_DEBUG 
LFLAGS.default.debug.	= $(LFLAGS.default.debug)
LFLAGS.default.ferret.	= $(LFLAGS.default.ferret)
LFLAGS.default.profile.	= $(LFLAGS.default.profile) -g -pg -O
LFLAGS.default.release.EXECUTABLE = -s 
LFLAGS.default.ferret.EXECUTABLE =  
LFLAGS.default.release.	= $(LFLAGS.default.release) -O \
			$(LFLAGS.default.release.$(BUILD_TYPE))
#ifeq (1,$(MAP))
LFLAGS += -Wl,-Map,$(TARGET_SPEC.EXECUTABLE).map
#endif

# support memory modeling
ifeq (1,$(HAVE_XTENSA))
ifeq (1,$(MEM_MODEL))
LFLAGS += -mlsp=sim-local
CFLAGS += -DMEM_MODEL=1
else
ifeq (2,$(MEM_MODEL))
LFLAGS += -mlsp=sim
CFLAGS += -DMEM_MODEL=2
else
endif
endif
endif

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
CFLAGSAVETEMP.1		= $(firstword $(CFLAGSAVETEMP.$(PLATFORM)) -save-temps)

LFLAGS_DEPENDENT	+= $(LFLAGS.$(PLATFORM))
LFLAGS.debug		= $(LFLAGS.debug.$(PLATFORM)) \
			$(LFLAGS.default.debug.$(LFLAGS.debug.$(PLATFORM)))
LFLAGS.profile		= $(LFLAGS.profile.$(PLATFORM)) \
			$(LFLAGS.default.profile.$(LFLAGS.profile.$(PLATFORM)))
LFLAGS.release		= $(LFLAGS.release.$(PLATFORM)) \
			$(LFLAGS.default.release.$(LFLAGS.release.$(PLATFORM)))
LFLAGS.ferret		= $(LFLAGS.ferret.$(PLATFORM)) \
			$(LFLAGS.default.ferret.$(LFLAGS.ferret.$(PLATFORM))) -Wl,-u -Wl,malloc -lferret
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
CFLAGS_DEPENDENT += $(CFLAGS.$(CONFIG))
LFLAGS_DEPENDENT += $(LFLAGS.$(CONFIG)) $(LFLAGS.$(BUILD_TYPE).$(SOLID_LIBRARY))
ifeq (0, $(CPP))
CFLAGS_DEPENDENT += -std=c99
endif
CFLAGS_DEPENDENT += $(IDIRS_OPT) -c -Dasm=__asm__
CFLAGS_DEPENDENT += -Wall -Wno-error 
#CFLAGS_DEPENDENT += -Wall -Werror 
CFLAGS_DEPENDENT += '-DFULL_VECTOR_DIR="./../../vectors_full"'
CFLAGS_DEPENDENT += '-DBRIEF_VECTOR_DIR="./../../vectors_brief"'

###############################################################
##                     Create target spec                    ##
###############################################################

#TARGET_SUFFIX.win32	= .exe

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

$(foreach lib,$(ALL_DLIBS),$(lib)): $(EXTERNAL_PROJECTS)

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

$(INTERMEDIATE_DIR)/$(TARGET).ld: $(LD_LIST.$(SOLID_LIBRARY))
	$(shell $(RM) $@) $(foreach mod,$(LD_LIST.$(SOLID_LIBRARY)),$(shell $(ECHO) $(mod)>>$@))

$(TARGET_SPEC.LIBRARY): $(EXTERNAL_PROJECTS) $(LD_LIST.$(SOLID_LIBRARY)) $(INTERMEDIATE_DIR)/$(TARGET).ld
	$(RM) $@
	$(AR) rcs $@ @$(INTERMEDIATE_DIR)/$(TARGET).ld

$(TARGET_SPEC.EXECUTABLE): $(EXTERNAL_PROJECTS) $(ALL_OBJECTS) $(ALL_DLIBS)
	$(LD) $(LFLAGS_DEPENDENT) $(LFLAGS) -o \
	$(TARGET_SPEC.EXECUTABLE) $(ALL_OBJECTS) $(ALL_LIBS)

all:    $(TARGET_SPEC.$(BUILD_TYPE))

$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)):
	+$(MAKE) -C $@ CONFIG=$(CONFIG) all

###############################################################
##                         Clean rules                       ##
###############################################################
$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/clean):
	$(MAKE) -C $(dir $@) CONFIG=$(CONFIG) clean
$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/cleanest):
	$(MAKE) -C $(dir $@) CONFIG=$(CONFIG) cleanest

clean:	$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/clean)
#clean:	
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.o)
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.$(BUILD_TYPE)))

cleanest:	$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/cleanest)
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.d)
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.o)
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.$(BUILD_TYPE)))

cleanlib:
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.LIBRARY)) 2>$(NULLDEV)
