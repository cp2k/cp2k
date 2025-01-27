PROJHOME := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
CP2KHOME := $(abspath $(PROJHOME)/../..)
OMP_TRACE := $(wildcard $(PROJHOME)/../base/openmp_trace.c)
ALL_HEADERS := $(shell find . -name "*.h") $(shell find ../offload/ -name "*.h")
ALL_OBJECTS :=  $(OMP_TRACE:.c=.o) \
        ../offload/offload_library.o \
        dbm_distribution.o \
        dbm_library.o \
        dbm_matrix.o \
        dbm_mempool.o \
        dbm_mpi.o \
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

CFLAGS := -g $(OPTFLAG) -march=native -Wall -Wextra -Wcast-qual

# Intel Compiler
ICX := $(shell which icx 2>/dev/null)
INTEL ?= $(if $(ICX),1,0)

ifeq ($(INTEL),0)
MKL_FCRTL := gf
LIBS += $(if $(OMPRT),-l$(OMPRT),-fopenmp)
CFLAGS += -fopenmp
else
MKL_FCRTL := intel
LIBS += $(if $(OMPRT),-l$(OMPRT),-qopenmp)
CFLAGS += -qopenmp
CC := icx
endif

# Some additional flags when relying on default-CC
ifeq ($(CC),)
CFLAGS += -Wno-vla-parameter
endif

NVCC := $(shell which nvcc 2>/dev/null)
NVARCH := sm_70
NVFLAGS := -g $(OPTFLAG) -lineinfo -arch $(NVARCH) -Wno-deprecated-gpu-targets -Xcompiler "$(CFLAGS)" -D__OFFLOAD_CUDA

HIPCC := $(shell which hipcc 2>/dev/null)
HIPARCH := gfx90a
HIPFLAGS := -fPIE -g $(OPTFLAG) --offload-arch=$(HIPARCH) -Wall -Wextra -Werror -I${ROCM_PATH}/include -D__OFFLOAD_HIP -D__HIP_PLATFORM_AMD__ -D__HIP_PLATFORM_AMD__

# prefer NVCC or HIPCC over OpenCL
ifneq ($(NVCC)$(HIPCC),)
NVCC_PATH := $(if $(NVCC),$(wildcard $(dir $(NVCC))/..))
CUDA_FILE := $(wildcard $(NVCC_PATH)/../cuda/include/cuda.h)
CUDA_PATH := $(if $(CUDA_FILE),$(NVCC_PATH)/../cuda,$(NVCC_PATH))
CUDA_LIBS := $(if $(wildcard $(CUDA_PATH)/lib64),lib64,lib)
OPENCL := 0
endif
ifneq ($(OPENCL),0)
OPENCL_OFFLOAD := $(firstword $(wildcard $(CP2KHOME)/lib/*/*/exts/dbcsr/libdbcsr.a))
OPENCL_GENKRNL := $(PROJHOME)/dbm_multiply_opencl.cl.h
ifneq ($(ICX),)
LIBS += -lifcoremt -limf
endif
endif

# Make foundational runtimes available
LIBS += -ldl -lstdc++ -lc -lm


# Make BLAS/LAPACK available
ifneq ($(MKLROOT),)
LIBS += -L$(MKLROOT)/lib/intel64 \
        -Wl,--start-group \
        -lmkl_$(MKL_FCRTL)_lp64 \
        -lmkl_core \
        -lmkl_sequential \
        -Wl,--end-group
else
LIBS += -lblas
endif

.PHONY : all clean

all: dbm_miniapp.x

clean:
	rm -fv $(ALL_OBJECTS) $(OPENCL_GENKRNL)

realclean: clean
	rm -fv dbm_miniapp.x


# Enable OpenCL when DBCSR library was prebuilt (assume __DBCSR_ACC).
ifneq ($(OPENCL_OFFLOAD),)
OPENCL_SRC := $(PROJHOME)/dbm_multiply_opencl.cl
OPENCL_CMN := $(wildcard $(CP2KHOME)/exts/dbcsr/src/acc/opencl/common/*.h)
OPENCL_GEN := $(CP2KHOME)/exts/dbcsr/src/acc/opencl/acc_opencl.sh
ALL_HEADERS += $(OPENCL_GENKRNL)

ALL_OBJECTS += dbm_multiply_gpu.o dbm_multiply_opencl.o
CFLAGS += -I$(CP2KHOME)/exts/dbcsr -D__OFFLOAD_OPENCL -D__DBCSR_ACC
LIBS += $(OPENCL_OFFLOAD) -lgfortran

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
ifneq (,$(wildcard $(OPENCL_ROOT)/include/sycl/CL/cl.h))
LIBS += -L$(OPENCL_ROOT)/$(if $(wildcard $(OPENCL_ROOT)/lib64),lib64,lib)
LIBS += -L$(OPENCL_ROOT)/compiler/lib/intel64 -lintlc
ifeq (,$(wildcard $(OPENCL_INC)))
CFLAGS += -I$(OPENCL_ROOT)/include/sycl
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

$(OPENCL_GENKRNL): $(OPENCL_GEN) $(OPENCL_SRC) $(OPENCL_CMN)
	$(OPENCL_GEN) -b 6 -p "" $(OPENCL_SRC) $@


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


# Make LIBXSMM available
ifneq ($(LIBXSMMROOT),)
CFLAGS += -D__LIBXSMM -I$(LIBXSMMROOT)/include
LIBS += \
        $(LIBXSMMROOT)/lib/libxsmm.a $(LIBXSMMROOT)/lib/libxsmmext.a \
        -lpthread -lrt
endif


# Build with MPI when mpicc compiler is present.
ifneq ($(MPICC),)
%.o: %.c $(OMP_TRACE) $(ALL_HEADERS)
	cd $(dir $<); $(MPICC) -c -std=c11 $(CFLAGS) -D__parallel $(notdir $<)

dbm_miniapp.x: dbm_miniapp.o $(ALL_OBJECTS)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LIBS)


# Build without MPI otherwise.
else
%.o: %.c $(OMP_TRACE) $(ALL_HEADERS)
	cd $(dir $<); $(CC) -c -std=c11 $(CFLAGS) $(notdir $<)

dbm_miniapp.x: dbm_miniapp.o $(ALL_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

endif
#EOF
