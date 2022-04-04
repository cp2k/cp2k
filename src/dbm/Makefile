.PHONY : all clean

all: dbm_miniapp.x

clean:
	rm -fv *.o */*.o ../offload/*.o

realclean: clean
	rm -fv *.x

CFLAGS := -fopenmp -g -O3 -march=native -Wall -Wextra
NVARCH := sm_70
NVFLAGS := -g -O3 -lineinfo -arch $(NVARCH) -Wno-deprecated-gpu-targets -Xcompiler "$(CFLAGS)" -D__OFFLOAD_CUDA

ifneq ($(MKLROOT),)
LIBS += \
        -Wl,--start-group \
        $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a \
        $(MKLROOT)/lib/intel64/libmkl_core.a \
        $(MKLROOT)/lib/intel64/libmkl_sequential.a \
        -Wl,--end-group -ldl
else
LIBS += -lblas
endif

LIBS += -lstdc++ -lm

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
NVCC := $(shell which nvcc)
ifneq ($(NVCC),)
LIBS += -lcudart -lcuda -lcublas -L${CUDA_PATH}/lib64
CFLAGS += -I${CUDA_PATH}/include -D__OFFLOAD_CUDA
ALL_OBJECTS += dbm_multiply_gpu.o

%.o: %.cu $(ALL_HEADERS)
	cd $(dir $<); $(NVCC) -c $(NVFLAGS) $(notdir $<)
endif

%.o: %.c $(ALL_HEADERS)
	cd $(dir $<); $(CC) -c -std=c11 $(CFLAGS) $(notdir $<)

dbm_miniapp.x: dbm_miniapp.o $(ALL_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

#EOF
