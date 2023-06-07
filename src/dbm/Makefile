.PHONY : all clean

all: dbm_miniapp.x

clean:
	rm -fv *.o */*.o ../offload/*.o

realclean: clean
	rm -fv *.x

CFLAGS := -fopenmp -g -O3 -march=native -Wall -Wextra -Wno-vla-parameter

NVCC := $(shell which nvcc)
NVARCH := sm_70
NVFLAGS := -g -O3 -lineinfo -arch $(NVARCH) -Wno-deprecated-gpu-targets -Xcompiler "$(CFLAGS)" -D__OFFLOAD_CUDA

HIPCC := $(shell which hipcc)
HIPARCH := gfx90a
HIPFLAGS := -fPIE -g -O3 --offload-arch=$(HIPARCH) -Wall -Wextra -Werror -I${ROCM_PATH}/include -D__OFFLOAD_HIP -D__HIP_PLATFORM_AMD__ -D__HIP_PLATFORM_AMD__

ifneq ($(LIBXSMMROOT),)
LIBS += \
        $(LIBXSMMROOT)/lib/libxsmm.a $(LIBXSMMROOT)/lib/libxsmmext.a \
        -lpthread -lrt
CFLAGS += -D__LIBXSMM -I$(LIBXSMMROOT)/include
endif

ifneq ($(MKLROOT),)
LIBS += \
        -L$(MKLROOT)/lib/intel64 \
        -Wl,--start-group \
        -lmkl_gf_lp64 \
        -lmkl_core \
        -lmkl_sequential \
        -Wl,--end-group
else
LIBS += -lblas
endif

LIBS += -fopenmp -ldl -lstdc++ -lc -lm

ALL_HEADERS := $(shell find . -name "*.h") $(shell find ../offload/ -name "*.h")
ALL_OBJECTS := ../offload/offload_library.o \
        dbm_distribution.o \
        dbm_library.o \
        dbm_matrix.o \
        dbm_mempool.o \
        dbm_mpi.o \
        dbm_multiply.o \
        dbm_multiply_comm.o \
        dbm_multiply_cpu.o \
        dbm_shard.o


# Enable Cuda when nvcc compiler is present.
ifneq ($(NVCC),)
LIBS += -lcudart -lcuda -lcublas -L${CUDA_PATH}/lib64
CFLAGS += -I${CUDA_PATH}/include -D__OFFLOAD_CUDA
ALL_OBJECTS += dbm_multiply_gpu.o dbm_multiply_gpu_kernel.o

%.o: %.cu $(ALL_HEADERS)
	cd $(dir $<); $(NVCC) -c $(NVFLAGS) $(notdir $<)


# Enable HIP/ROCm when hipcc compiler is present.
else ifneq ($(HIPCC),)
LIBS += -L${ROCM_PATH}/lib -lamdhip64 -lhipfft -lhipblas -lrocblas
CFLAGS += -I${ROCM_PATH}/include -D__OFFLOAD_HIP -D__HIP_PLATFORM_AMD__
ALL_OBJECTS += dbm_multiply_gpu.o

%.o: %.cu $(ALL_HEADERS)
	cd $(dir $<); $(HIPCC) -c $(HIPFLAGS) $(notdir $<)
endif


# Build with MPI when mpicc compiler is present.
ifneq ($(MPICC),)
%.o: %.c $(ALL_HEADERS)
	cd $(dir $<); $(MPICC) -c -std=c11 $(CFLAGS) -D__parallel $(notdir $<)

dbm_miniapp.x: dbm_miniapp.o $(ALL_OBJECTS)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LIBS)


# Build without MPI otherwise.
else
%.o: %.c $(ALL_HEADERS)
	cd $(dir $<); $(CC) -c -std=c11 $(CFLAGS) $(notdir $<)

dbm_miniapp.x: dbm_miniapp.o $(ALL_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

endif
#EOF
