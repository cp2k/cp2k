.PHONY : all clean

all: grid_miniapp.x grid_unittest.x

clean:
	rm -fv *.o */*.o *.x ../offload/*.o

CFLAGS := -fopenmp -g -O3 -march=native -Wall -Wextra -Wno-vla-parameter
NVFLAGS := -g -O3 -lineinfo -arch sm_70 -Wno-deprecated-gpu-targets -Xcompiler "$(CFLAGS)" -D__OFFLOAD_CUDA
LIBS := -lm -lblas

ALL_HEADERS := $(shell find . -name "*.h") $(shell find ../offload/ -name "*.h")
ALL_OBJECTS := ../offload/offload_buffer.o \
        ../offload/offload_library.o \
        grid_replay.o \
        grid_task_list.o \
        grid_collocate.o \
        grid_integrate.o \
        grid_multigrid.o \
        common/grid_library.o \
        common/grid_mpi.o \
        common/grid_basis_set.o \
        common/grid_sphere_cache.o \
        ref/grid_ref_task_list.o \
        ref/grid_ref_collocate.o \
        ref/grid_ref_integrate.o \
        ref/grid_ref_prepare_pab.o \
        cpu/grid_cpu_task_list.o \
        cpu/grid_cpu_collocate.o \
        cpu/grid_cpu_integrate.o \
        cpu/grid_cpu_prepare_pab.o \
        dgemm/grid_dgemm_context.o \
        dgemm/grid_dgemm_coefficients.o \
        dgemm/grid_dgemm_collocate.o \
        dgemm/grid_dgemm_integrate.o \
        dgemm/grid_dgemm_non_orthorombic_corrections.o \
        dgemm/grid_dgemm_utils.o \
        dgemm/grid_dgemm_collocation_integration.o \
        dgemm/grid_dgemm_context.o \
        dgemm/grid_dgemm_prepare_pab.o \
        dgemm/grid_dgemm_tensor_local.o

# Enable Cuda when nvcc compiler is present.
NVCC := $(shell which nvcc)
ifneq ($(NVCC),)
LIBS += -lcudart -lcuda -lcublas -L${CUDA_PATH}/lib64
CFLAGS += -I${CUDA_PATH}/include -D__OFFLOAD_CUDA
ALL_OBJECTS += gpu/grid_gpu_task_list.o \
        gpu/grid_gpu_collocate.o \
        gpu/grid_gpu_integrate.o

%.o: %.cu $(ALL_HEADERS)
	cd $(dir $<); $(NVCC) -c $(NVFLAGS) $(notdir $<)
endif

%.o: %.c $(ALL_HEADERS)
	cd $(dir $<); $(CC) -c -std=c11 $(CFLAGS) $(notdir $<)

grid_miniapp.x: grid_miniapp.o $(ALL_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

grid_unittest.x: grid_unittest.o $(ALL_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

#EOF
