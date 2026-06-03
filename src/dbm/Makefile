PROJHOME := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
CP2KHOME := $(abspath $(PROJHOME)/../..)
OMP_TRACE := $(wildcard $(PROJHOME)/../base/openmp_trace.c)
ALL_HEADERS := $(shell find . -name "*.h") $(shell find ../offload/ -name "*.h")
ALL_OBJECTS :=  $(OMP_TRACE:.c=.o) \
        ../offload/offload_library.o \
        ../offload/offload_mempool.o \
        ../mpiwrap/cp_mpi.o \
        dbm_distribution.o \
        dbm_library.o \
        dbm_matrix.o \
        dbm_multiply.o \
        dbm_multiply_comm.o \
        dbm_multiply_cpu.o \
        dbm_shard.o

# Optimization level
ifeq (,$(filter-out 0,$(DBG)))
OPT ?= 3
else
OPT ?= 0
endif

# Optimization flag derived from OPT flag
OPTFLAG ?= -O$(patsubst O%,%,$(OPT))

# Common compiler flags
CFLAGS := -fPIC -g $(OPTFLAG) -march=native -Wall -Wextra

PEDANTIC ?= 0
ifneq ($(PEDANTIC),0)
CFLAGS += -Wcast-qual -pedantic -Wno-overlength-strings
endif

USM ?= 0
ifneq ($(USM),0)
CFLAGS += -D__OFFLOAD_UNIFIED_MEMORY
endif

# Intel Compiler
ICX := $(shell which icx 2>/dev/null)
INTEL ?= $(if $(ICX),$(if $(filter-out 0,$(GNU)),0,1),0)

# Build with MPI when MPICC is given
MPI ?= $(if $(MPICC),1,0)

ifeq ($(INTEL),0)
MKL_FCRTL := gf
LIBS += $(if $(OMPRT),-l$(OMPRT),-fopenmp)
CFLAGS += -fopenmp -Wno-vla-parameter
CC := $(if $(filter-out 0,$(MPI)),mpicc,gcc)
else
MKL_FCRTL := intel
LIBS += $(if $(OMPRT),-fopenmp -Wno-recommended-option,-qopenmp)
CFLAGS += $(if $(OMPRT),-fopenmp -Wno-recommended-option,-qopenmp)
CC := $(if $(filter-out 0,$(MPI)),mpiicx,icx)
endif

ifneq ($(MPI),0)
CFLAGS += -D__parallel
ifneq ($(MPICC),)
CC := $(MPICC)
endif
endif

# Make LIBXS available (auto-detect as sibling of CP2KHOME).
ifneq ($(LIBXS),0)
LIBXSROOT ?= $(wildcard $(CP2KHOME)/../libxs)
ifneq ($(LIBXSROOT),)
CFLAGS += -D__LIBXS -D__BLAS -I$(LIBXSROOT)/include
LIBS += $(LIBXSROOT)/lib/libxs.a
# LIBXSMM headers enable JIT kernel dispatch inside LIBXS.
LIBXSMMROOT ?= $(wildcard $(LIBXSROOT)/../libxsmm)
ifneq ($(LIBXSMMROOT),)
ifneq ($(wildcard $(LIBXSMMROOT)/lib/libxsmm.a),)
CFLAGS += -D__LIBXSMM -I$(LIBXSMMROOT)/include
LIBS += $(LIBXSMMROOT)/lib/libxsmm.a
endif
endif
else
OPENCL := 0
endif
endif

# Make LIBXSTREAM available (auto-detect as sibling of CP2KHOME).
ifneq ($(LIBXSTREAM),0)
ifneq ($(OPENCL),0)
ifneq ($(LIBXS),0)
LIBXSTREAMROOT ?= $(wildcard $(CP2KHOME)/../libxstream)
ifneq ($(LIBXSTREAMROOT),)
CFLAGS += -D__LIBXSTREAM -I$(LIBXSTREAMROOT)/include
LIBS += $(LIBXSTREAMROOT)/lib/libxstream.a
else
OPENCL := 0
endif
endif
endif
endif

ifeq ($(filter-out 0,$(OPENCL)),)
NVCC := $(shell which nvcc 2>/dev/null)
NVARCH := sm_70
NVARCHX := $(if $(NVARCH),-arch $(NVARCH))
NVFLAGS := -g $(OPTFLAG) -lineinfo $(NVARCHX) -Wno-deprecated-gpu-targets -Xcompiler "$(CFLAGS)" -D__OFFLOAD_CUDA

HIPCC := $(shell which hipcc 2>/dev/null)
HIPARCH := gfx90a
ROCM_PATH_ENV := $(shell echo "$${ROCM_PATH}")
ROCM_PATH := $(if $(ROCM_PATH_ENV),$(ROCM_PATH_ENV),/opt/rocm)
HIPFLAGS := -fPIE -g $(OPTFLAG) --offload-arch=$(HIPARCH) -Wall -Wextra -Werror -I${ROCM_PATH}/include -D__OFFLOAD_HIP -D__HIP_PLATFORM_AMD__
endif

# prefer NVCC or HIPCC over OpenCL
ifneq ($(NVCC)$(HIPCC),)
NVCC_PATH := $(if $(NVCC),$(wildcard $(dir $(NVCC))/..))
CUDA_FILE := $(wildcard $(NVCC_PATH)/../cuda/include/cuda.h)
CUDA_PATH := $(if $(CUDA_FILE),$(NVCC_PATH)/../cuda,$(NVCC_PATH))
CUDA_LIBS := $(if $(wildcard $(CUDA_PATH)/lib64),lib64,lib)
OPENCL := 0
endif
ifneq ($(OPENCL),0)
OPENCL_GENKRNL := $(PROJHOME)/dbm_multiply_opencl.cl.h
ifneq ($(ICX),)
LIBS += -lifcoremt -limf
endif
endif

# Make BLAS/LAPACK available
ifneq ($(MKLROOT),)
LIBS += -L$(MKLROOT)/lib/intel64 \
        -Wl,--start-group \
        -lmkl_$(MKL_FCRTL)_lp64 \
        -lmkl_core \
        -lmkl_sequential \
        -Wl,--end-group
ifneq ($(MKL),0)
CFLAGS += -D__MKL -I$(MKLROOT)/include
endif
else ifneq ($(OPENBLAS),0)
LIBS += -lopenblas
else
LIBS += -lblas
endif

# Make foundational runtimes available
LIBS += -ldl -lstdc++ -lc -lm

.PHONY : all clean

all: dbm_miniapp.x

clean:
	rm -fv $(ALL_OBJECTS) $(OPENCL_GENKRNL) \
        dbm_multiply_gpu.o \
        dbm_multiply_gpu_kernel.o \
        dbm_multiply_opencl.o \
        dbm_miniapp.o

realclean: clean
	rm -fv dbm_miniapp.x


# Enable OpenCL offload backend.
ifneq ($(OPENCL),0)
OPENCL_SRC := $(PROJHOME)/dbm_multiply_opencl.cl
OPENCL_CMN := $(wildcard $(LIBXSTREAMROOT)/include/opencl/*.h)
OPENCL_GEN := $(LIBXSTREAMROOT)/scripts/tool_opencl.sh
OPENCL_ARG := -b 6 -p "" -I $(LIBXSTREAMROOT)/include

ifneq (,$(filter-out 0,$(DBG)))
OPENCL_ARG += -d
endif

ALL_HEADERS += $(OPENCL_GENKRNL)
ALL_OBJECTS += dbm_multiply_gpu.o dbm_multiply_opencl.o
CFLAGS += -D__OFFLOAD_OPENCL -I$(CP2KHOME)/exts/dbcsr
LIBS += -lgfortran

ifneq (Darwin,$(shell uname))
OPENCL_LIB := $(shell ldconfig -p 2>/dev/null | grep -m1 OpenCL | rev | cut -d' ' -f1 | rev)
ifeq (,$(OPENCL_LIB))
OPENCL_LIB := /usr/lib/x86_64-linux-gnu/libOpenCL.so.1
endif
ifneq (,$(CUDA_PATH))
ifeq (,$(wildcard $(OPENCL_INC)))
CFLAGS += -I$(CUDA_PATH)/include
endif
ifeq (,$(wildcard $(OPENCL_LIB)))
LIBS += -L$(CUDA_PATH)/$(CUDA_LIBS)
LIBS += -Wl,-rpath=$(CUDA_PATH)/$(CUDA_LIBS)
endif
else ifneq (,$(wildcard $(OPENCL_ROOT)/include/CL/cl.h))
ifeq (,$(wildcard $(OPENCL_INC)))
CFLAGS += -I$(OPENCL_ROOT)/include
endif
LIBS += -L$(OPENCL_ROOT)/$(if $(wildcard $(OPENCL_ROOT)/lib64),lib64,lib)
else ifneq (,$(ICX))
OPENCL_ROOT := $(abspath $(dir $(ICX))/..)
CLINC := $(wildcard $(OPENCL_ROOT)/include/sycl/CL/cl.h $(OPENCL_ROOT)/include/CL/cl.h)
ifneq (,$(CLINC))
LIBS += -L$(OPENCL_ROOT)/$(if $(wildcard $(OPENCL_ROOT)/lib64),lib64,lib)
LIBS += -L$(OPENCL_ROOT)/compiler/lib/intel64 -lintlc
ifeq (,$(wildcard $(OPENCL_INC)))
CFLAGS += -I$(abspath $(dir $(firstword $(CLINC)))/..)
endif
endif
endif
# OPENCL_INC: directory containing CL/cl.h.
ifneq (,$(wildcard $(OPENCL_INC)))
CFLAGS += -I$(OPENCL_INC)
endif
# OPENCL_LIB: file/library to be linked
ifneq (,$(wildcard $(OPENCL_LIB)))
LIBS += $(OPENCL_LIB)
else
LIBS += -l:libOpenCL.so.1
endif
else # macOS
LIBS += -framework OpenCL
endif

dbm_miniapp.x: $(OPENCL_OFFLOAD)
$(OPENCL_GENKRNL): $(OPENCL_GEN) $(OPENCL_SRC) $(OPENCL_CMN)
	$(OPENCL_GEN) $(OPENCL_ARG) $(OPENCL_SRC) $@


# Enable Cuda when nvcc compiler is present.
else ifneq ($(NVCC),)
ALL_OBJECTS += dbm_multiply_gpu.o dbm_multiply_gpu_kernel.o
CFLAGS += -I${CUDA_PATH}/include -D__OFFLOAD_CUDA
LIBS += -lcudart -lcuda -lcublas -L${CUDA_PATH}/$(CUDA_LIBS)
ifneq ($(wildcard $(NVCC_PATH)/../math_libs/$(CUDA_LIBS)),)
LIBS += -L$(NVCC_PATH)/../math_libs/$(CUDA_LIBS)
endif

%.o: %.cu $(ALL_HEADERS)
	cd $(dir $<); $(NVCC) -c $(NVFLAGS) $(notdir $<)


# Enable HIP/ROCm when hipcc compiler is present.
else ifneq ($(HIPCC),)
ALL_OBJECTS += dbm_multiply_gpu.o dbm_multiply_gpu_kernel.o
CFLAGS += -I${ROCM_PATH}/include -D__OFFLOAD_HIP -D__HIP_PLATFORM_AMD__
LIBS += -L${ROCM_PATH}/lib -lamdhip64 -lhipfft -lhipblas -lrocblas

%.o: %.cu $(ALL_HEADERS)
	cd $(dir $<); $(HIPCC) -c $(HIPFLAGS) $(notdir $<)
endif

# LIBXS must be linked last (depended upon by LIBXSTREAM)
ifneq ($(LIBXSROOT),)
LIBS += $(LIBXSROOT)/lib/libxs.a
endif

%.o: %.c $(OMP_TRACE) $(ALL_HEADERS)
	cd $(dir $<); $(CC) -c -std=c11 $(CFLAGS) $(notdir $<)

dbm_miniapp.x: dbm_miniapp.o $(ALL_OBJECTS)
	$(CC) -o $@ $^ $(LIBS)

#EOF
