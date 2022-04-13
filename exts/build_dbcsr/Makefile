# This Makefile is only used for the inclusion of the DBCSR library into CP2K.
# It is not supported by DBCSR development anymore.
############################################################
#### DO NOT CHANGE AS PART OF THE DBCSR DEVELOPMENT!!!! ####
############################################################
# CP2K team will update it, following the compilation flags suggested in the
# CP2K installation file.
# For this reason, this Makefile must be compatible with the CP2K
# compilation and tested within CP2K.
#
# The ARCH file can set the variables:
#
# CC  => C compiler, e.g. gcc or mpicc
# FC  => Fortran compiler, e.g. gfortran or mpifort
# LD  => Linker, e.g. gfortran or mpifort
# AR  => Archive command, e.g. ar -r
# CXXFLAGS  => C++ compilation flags
# CFLAGS    => C compilation flags
# FCFLAGS   => Fortran compilation flags
# LDFLAGS   => Linker flags
# LIBS      => Libraries
# ACC       => ACC can be nvcc (CUDA), mpicxx/mpicc or gcc/g++ (HIP), or mpicc/gcc (OpenCL)
# ACCFLAGS  => ACC flags
# USE_ACCEL => hip, cuda, or opencl
# GPUVER    =>
#          - for CUDA, possible values correspond to NVIDIA GPUs:
#            possible values are K20X, K40, K80, P100, V100
#          - for HIP, possible values correspond to NVIDIA and AMD GPUs:
#            possible values are K20X, K40, K80, P100, V100, Mi50, Mi100
#          - for OpenCL, GPUVER maps to a file with tuned parameters:
#            src/acc/opencl/smm/params/tune_multiply_GPUVER.csv, or
#            src/acc/opencl/smm/tune_multiply.csv (default)
#
# Libraries for accelerator:
#    - e.g. for CUDA:   LIBS += -lstdc++ -lcudart -lnvrtc -lcuda -lcublas
#    - e.g. for HIP:    LIBS += -lstdc++ -lhiprtc -lhipblas
#    - e.g. for OpenCL: LIBS += -lOpenCL

#
SHELL = /bin/sh
#
# the home dir is taken from the current directory
#
DBCSRHOME    := $(CURDIR)
DBCSRCP2K    := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
MAKEFILE     := $(DBCSRCP2K)/Makefile
LIBDIR       := $(DBCSRHOME)/lib
OBJDIR       := $(DBCSRHOME)/obj
TOOLSDIR     := $(DBCSRHOME)/tools
FYPPEXE      := $(TOOLSDIR)/build_utils/fypp/bin/fypp
SRCDIR       := $(DBCSRHOME)/src
TESTSDIR     := $(DBCSRHOME)/tests
ARCHFILE     :=

PYTHON       := /usr/bin/env python3

# Default Target ============================================================
LIBNAME      := dbcsr
LIBRARY      := lib$(LIBNAME)
default_target: $(LIBRARY)

# Read the version ==========================================================
include $(DBCSRHOME)/VERSION
ifeq ($(DATE),)
  DATE = "Development Version"
endif

# Read the configuration ====================================================
ifeq (,$(ARCHFILE))
$(error ARCH file must be provided: ARCHFILE=<CP2K ARCH file>)
endif
ifeq (, $(wildcard $(ARCHFILE)))
  $(error ARCH file `$(ARCHFILE)` not found)
endif
include $(ARCHFILE)

# Set the ARCH version ======================================================
ifeq ($(GPUVER),K20X)
  ARCH_NUMBER = 35
else ifeq ($(GPUVER),K40)
  ARCH_NUMBER = 35
else ifeq ($(GPUVER),K80)
  ARCH_NUMBER = 37
else ifeq ($(GPUVER),P100)
  ARCH_NUMBER = 60
else ifeq ($(GPUVER),V100)
  ARCH_NUMBER = 70
else ifeq ($(GPUVER),A100)
  # TODO: update for A100 tuned parameters
  override GPUVER := V100
  ARCH_NUMBER = 80
else ifeq ($(GPUVER),Mi50)
  ARCH_NUMBER = gfx906
else ifeq ($(GPUVER),Mi100)
  ARCH_NUMBER = gfx908
else ifeq (,$(ARCH_NUMBER))
ifneq (opencl,$(USE_ACCEL))
  $(error Unknown ARCH_NUMBER since GPUVER="$(GPUVER)" is not recognized)
endif
endif

# Set ACCFLAGS ==============================================================

# If compiling with nvcc
ifeq (cuda,$(USE_ACCEL))
  override ACCFLAGS += -D__CUDA
  FCFLAGS += -D__CUDA
  CXXFLAGS += -D__CUDA
  #if "-arch" has not yet been set in ACCFLAGS
  ifeq (,$(findstring -arch,$(ACCFLAGS)))
    override ACCFLAGS += -arch sm_$(ARCH_NUMBER)
  endif
  ifeq (,$(findstring -Xcompiler,$(ACCFLAGS)))
    override ACCFLAGS += -Xcompiler="$(CXXFLAGS)"
  endif
# If compiling with hipcc
else ifeq (hip,$(USE_ACCEL))
  override ACCFLAGS += -D__HIP
  FCFLAGS += -D__HIP
# OpenCL backend
else ifeq (opencl,$(USE_ACCEL))
  CFLAGS = $(ACCFLAGS) -D__OPENCL
  # LIBXSMM is eventually included privately/quietly (header-only)
  LIBXSMMROOT := $(if $(LIBXSMMROOT),$(LIBXSMMROOT),$(if $(LIBXSMM_INC),$(wildcard $(LIBXSMM_INC)/..)))
  LIBXSMMROOT := $(if $(LIBXSMMROOT),$(LIBXSMMROOT),$(wildcard ../../../libxsmm))
  LIBXSMMROOT := $(if $(LIBXSMMROOT),$(LIBXSMMROOT),$(wildcard $(HOME)/libxsmm))
  ifneq (,$(LIBXSMMROOT))
    CFLAGS += -I$(LIBXSMMROOT)/include
  endif
endif

# Set the configuration =====================================================
ifneq ($(LD_SHARED),)
  ARCHIVE_EXT := .so
else
  ARCHIVE_EXT := .a
endif

# Declare PHONY targets =====================================================
.PHONY : dirs makedep \
         default_target $(LIBRARY) \
         clean version

# Discover files and directories ============================================
ALL_SRC_DIRS       := $(shell find $(SRCDIR) -type d | awk '{printf("%s:",$$1)}')
LIBSMM_ACC_DIR     := $(shell cd $(SRCDIR) ; find . -type d -name "libsmm_acc")
LIBSMM_ACC_ABS_DIR := $(shell find $(SRCDIR) -type d -name "libsmm_acc")

ALL_PKG_FILES := $(shell find $(SRCDIR) -name "PACKAGE")
OBJ_SRC_FILES  = $(shell cd $(SRCDIR); find . ! -name "dbcsr_api_c.F" ! -name "dbcsr_tensor_api_c.F" -name "*.F")
OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . -type f -name "*.c")

# Script used to create dependency/file which depends on GPUVER
ACC_MAKEDEP := $(wildcard $(LIBSMM_ACC_ABS_DIR)/../acc_makedep.sh)

# if compiling with GPU acceleration
ifneq ($(ACC),)
ifeq (cuda,$(USE_ACCEL))
  # All *.cpp files belong to the accelerator backend common between CUDA/HIP
  OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . -type f ! -name "acc_cuda.cpp" ! -name "acc_hip.cpp" -name "*.cpp")
  OBJ_SRC_FILES += $(LIBSMM_ACC_DIR)/../cuda/acc_cuda.cpp
  # Exclude autotuning files
  OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . -type f ! -name "tune_*_exe*_part*.cu" ! -name "tune_*_exe*_main*.cu"  -name "*.cu")
else ifeq (hip,$(USE_ACCEL))
  # All *.cpp files belong to the accelerator backend common between CUDA/HIP
  OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . -type f ! -name "acc_cuda.cpp" ! -name "acc_hip.cpp" -name "*.cpp")
  OBJ_SRC_FILES += $(LIBSMM_ACC_DIR)/../hip/acc_hip.cpp
# OpenCL backend: OBJ_SRC_FILES already includes all *.c files
endif
endif

# Include also source files which won't compile into an object file
ALL_SRC_FILES  = $(strip $(subst $(NULL) .,$(NULL) $(SRCDIR),$(NULL) $(OBJ_SRC_FILES)))
ALL_SRC_FILES += $(shell find $(SRCDIR) -name "*.f90")
ALL_SRC_FILES += $(shell find $(SRCDIR) -name "*.h")
ALL_SRC_FILES += $(shell find $(SRCDIR) -name "*.hpp")

# stage 1: create dirs and run makedep.py.
#          Afterwards, call make recursively again with -C $(OBJDIR) and INCLUDE_DEPS=true
ifeq ($(INCLUDE_DEPS),)
$(LIBRARY): dirs makedep
	@+$(MAKE) --no-print-directory -C $(OBJDIR) -f $(MAKEFILE) $(LIBDIR)/$(LIBRARY)$(ARCHIVE_EXT) INCLUDE_DEPS=true DBCSRHOME=$(DBCSRHOME)
dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(LIBDIR)
version:
	@echo "DBCSR Version: "$(MAJOR)"."$(MINOR)"."$(PATCH)" ("$(DATE)")"
else
# stage 2: Include $(OBJDIR)/all.dep, expand target all, and get list of dependencies.

# Check if FYPP is available  ===============================================
ifeq (, $(shell which $(FYPPEXE) 2>/dev/null ))
  $(error No FYPP submodule available, please read README.md on how to properly download DBCSR)
endif

endif

clean:
	rm -f $(LIBSMM_ACC_ABS_DIR)/parameters.h $(LIBSMM_ACC_ABS_DIR)/smm_acc_kernels.h $(LIBSMM_ACC_ABS_DIR)/*.so
	rm -f $(LIBSMM_ACC_ABS_DIR)/../opencl/smm/opencl_kernels.h
	rm -f $(TESTSDIR)/libsmm_acc_unittest_multiply.cpp
	rm -f $(TESTSDIR)/libsmm_acc_timer_multiply.cpp
	rm -rf $(OBJDIR)

# Libsmm_acc stuff ==========================================================
ifneq ($(ACC),)
ifneq (opencl,$(USE_ACCEL))
ACC_KERNEL := $(wildcard $(LIBSMM_ACC_ABS_DIR)/kernels/*.h)
ACC_PARAMS := $(wildcard $(LIBSMM_ACC_ABS_DIR)/parameters/parameters_$(GPUVER).json)
ACC_PARDEP := $(if $(ACC_MAKEDEP),$(shell $(ACC_MAKEDEP) $(LIBSMM_ACC_ABS_DIR)/.with_gpu $(GPUVER)))
$(LIBSMM_ACC_ABS_DIR)/parameters.h: $(LIBSMM_ACC_ABS_DIR)/generate_parameters.py $(ACC_PARAMS) $(ACC_PARDEP)
	cd $(LIBSMM_ACC_ABS_DIR); $(PYTHON) generate_parameters.py --gpu_version=$(GPUVER)
$(LIBSMM_ACC_ABS_DIR)/smm_acc_kernels.h: $(LIBSMM_ACC_ABS_DIR)/generate_kernels.py $(ACC_KERNEL)
	cd $(LIBSMM_ACC_ABS_DIR); $(PYTHON) generate_kernels.py
endif
endif

# automatic dependency generation ===========================================
MODDEPS     = "lower"
MAKEDEPMODE = "normal"
ifeq ($(HACKDEP),yes)
MAKEDEPMODE = "hackdep"
endif

# this happens at stage 1
makedep: $(ALL_SRC_FILES) $(ALL_PKG_FILES) dirs
ifeq ($(LD_SHARED),)
	@echo "Removing stale archives ... "
	@$(PYTHON) $(DBCSRCP2K)/check_archives.py $(firstword $(AR)) $(SRCDIR) $(LIBDIR)
endif
	@echo "Resolving dependencies ... "
	@$(PYTHON) $(DBCSRCP2K)/makedep.py $(OBJDIR)/all.dep dbcsr $(MODDEPS) $(MAKEDEPMODE) $(ARCHIVE_EXT) $(SRCDIR) $(OBJ_SRC_FILES)

# at stage 2, load the rules generated by makedep.py
ifeq ($(INCLUDE_DEPS), true)
  include $(OBJDIR)/all.dep
endif


# ================= Stuff need for compiling (stage 2) ======================
# These rules are executed in a recursive call to make -C $(OBJDIR)
# The change of $(CURDIR) allows to find targets without abs paths and vpaths.


### Slave rules ###
vpath %.F     $(ALL_SRC_DIRS)
vpath %.h     $(ALL_SRC_DIRS)
vpath %.hpp   $(ALL_SRC_DIRS)
vpath %.f90   $(ALL_SRC_DIRS)
vpath %.cu    $(ALL_SRC_DIRS)
vpath %.c     $(ALL_SRC_DIRS)
vpath %.cpp   $(ALL_SRC_DIRS)

# $(FCLOGPIPE) can be used to store compiler output, e.g. warnings, for each F-file separately.
# This is used e.g. by the convention checker.

FYPPFLAGS ?= -n

%.o: %.F
	$(PYTHON) $(FYPPEXE) $(FYPPFLAGS) $< $*.F90
	$(FC) -c $(FCFLAGS) -D__SHORT_FILE__="\"$(notdir $<)\"" -I'$(dir $<)' -I'$(SRCDIR)' $*.F90 $(FCLOGPIPE)

%.mod: %.o
	@true

%.o: %.c
	$(CC) -c $(CFLAGS) $<

# Compile the CUDA/HIP files
ifneq ($(ACC),)
ifneq (opencl,$(USE_ACCEL))
%.o: %.cpp
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
libsmm_acc.o: libsmm_acc.cpp parameters.h smm_acc_kernels.h
	$(ACC) -c $(ACCFLAGS) -DARCH_NUMBER=$(ARCH_NUMBER) $<
libsmm_acc_benchmark.o: libsmm_acc_benchmark.cpp parameters.h
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
libsmm_acc_init.o: libsmm_acc_init.cpp libsmm_acc_init.h parameters.h
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
endif
endif

# if compiling CUDA backend
ifeq (cuda,$(USE_ACCEL))
ifneq ($(ACC),)
%.o: %.cpp
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
acc_cuda.o: acc_cuda.cpp acc_cuda.h
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
%.o: %.cu
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
# if compiling HIP backend
else ifeq (hip,$(USE_ACCEL))
%.o: %.cpp
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
acc_hip.o: acc_hip.cpp acc_hip.h
	$(ACC) -c $(ACCFLAGS) -I'$(SRCDIR)' $<
endif
# if compiling OpenCL backend
else ifeq (opencl,$(USE_ACCEL))
OPENCL_KRNLGEN := $(LIBSMM_ACC_ABS_DIR)/../opencl/acc_opencl.sh
OPENCL_KERNELS := $(wildcard $(LIBSMM_ACC_ABS_DIR)/../opencl/smm/kernels/*.cl)
OPENCL_DEFAULT := $(wildcard $(LIBSMM_ACC_ABS_DIR)/../opencl/smm/tune_multiply.csv)
OPENCL_WITHGPU := $(wildcard $(LIBSMM_ACC_ABS_DIR)/../opencl/smm/params/tune_multiply_$(GPUVER).csv)
OPENCL_PARAMS := $(if $(OPENCL_WITHGPU),$(OPENCL_WITHGPU),$(OPENCL_DEFAULT))
OPENCL_PARDEP := $(if $(ACC_MAKEDEP),$(shell $(ACC_MAKEDEP) $(LIBSMM_ACC_ABS_DIR)/../opencl/smm/.with_gpu $(GPUVER)))
$(LIBSMM_ACC_ABS_DIR)/../opencl/smm/opencl_kernels.h: $(OPENCL_KRNLGEN) $(OPENCL_KERNELS) $(OPENCL_PARAMS) $(OPENCL_PARDEP)
	$(OPENCL_KRNLGEN) $(OPENCL_KERNELS) $(OPENCL_PARAMS) $@
opencl_libsmm.o: opencl_libsmm.c $(LIBSMM_ACC_ABS_DIR)/../opencl/smm/opencl_kernels.h
ifeq (Darwin,$(shell uname))
  LDFLAGS += -framework OpenCL
else
  # OpenCL include directory (cl.h not installed per "opencl-headers" package)
  ifeq (,$(CUDATOOLKIT_HOME))
    CUDATOOLKIT_HOME := $(NVSDKCOMPUTE_ROOT)
  endif
  ifeq (,$(CUDATOOLKIT_HOME))
    NVCC := $(call which,nvcc)
    CUDATOOLKIT_HOME := $(if $(NVCC),$(abspath $(dir $(NVCC))/..))
  endif
  ifneq (,$(CUDATOOLKIT_HOME))
    CFLAGS += -I$(CUDATOOLKIT_HOME)/include
  endif
endif
endif

$(LIBDIR)/%:
ifneq ($(LD_SHARED),)
	@echo "Creating shared library $@"
	@$(LD_SHARED) $(LDFLAGS) -o $(@:.a=.so) $^ $(LIBS)
else
	@echo "Updating archive $@"
	@$(AR) $@ $?
endif

#EOF
