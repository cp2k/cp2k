.PHONY : all clean

all: grid_miniapp.x grid_unittest.x

clean:
	rm -fv *.o */*.o *.x

CFLAGS := -fopenmp -g -O3 -march=native -Wall -Wextra
NVFLAGS := -g -O3 -lineinfo -arch sm_50 -Wno-deprecated-gpu-targets -Xcompiler "$(CFLAGS)" -D__GRID_CUDA
LIBS := -lm -lblas

ALL_HEADERS := $(shell find . -name "*.h")
ALL_OBJECTS := grid_replay.o \
							 grid_task_list.o \
							 common/grid_library.o \
							 common/grid_buffer.o \
							 common/grid_basis_set.o \
							 common/grid_sphere_cache.o \
							 ref/grid_ref_task_list.o \
							 ref/grid_ref_collocate.o \
							 ref/grid_ref_integrate.o \
							 ref/grid_ref_prepare_pab.o \
							 cpu/grid_context_cpu.o \
							 cpu/coefficients.o \
							 cpu/grid_collocate_dgemm.o \
							 cpu/grid_integrate_dgemm.o \
							 cpu/non_orthorombic_corrections.o \
							 cpu/utils.o \
							 cpu/collocation_integration.o \
							 cpu/grid_context_cpu.o \
							 cpu/grid_prepare_pab_dgemm.o \
							 cpu/tensor_local.o

# Enable Cuda when nvcc compiler is present.
NVCC := $(shell which nvcc)
ifneq ($(NVCC),)
LIBS += -lcudart -lcuda -lcublas -L${CUDA_PATH}/lib64
CFLAGS += -I${CUDA_PATH}/include -D__GRID_CUDA
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
