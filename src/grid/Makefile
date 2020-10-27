.PHONY : all clean

all: grid_collocate_miniapp.x grid_collocate_unittest.x

clean:
	rm -fv grid_*.o */grid_*.o grid_*.x

CFLAGS := -fopenmp -g -O3 -std=c99 -march=native -Wall -Wextra -Wpedantic -I${CUDA_PATH}/include -D__GRID_CUDA
NVFLAGS := -g -O3 -arch sm_50 -Wno-deprecated-gpu-targets -Xcompiler "-fopenmp -Wall -Wextra" -D__GRID_CUDA

LIBS := -lm -lblas -lcudart -lcuda -L${CUDA_PATH}/lib64

ALL_HEADERS := $(shell find . -name "*.h")

ALL_OBJECTS := grid_collocate_replay.o \
               grid_task_list.o \
               common/grid_library.o \
               common/grid_basis_set.o \
               common/grid_sphere_cache.o \
               cpu/grid_cpu_task_list.o \
               ref/grid_ref_task_list.o \
               ref/grid_ref_collocate.o \
               ref/grid_ref_integrate.o \
               ref/grid_ref_prepare_pab.o \
               gpu/grid_gpu_task_list.o \
               gpu/grid_gpu_collocate.o

%.o: %.c $(ALL_HEADERS)
	cd $(dir $<); $(CC) -c $(CFLAGS) $(notdir $<)

%.o: %.cu $(ALL_HEADERS)
	cd $(dir $<); nvcc -c $(NVFLAGS) $(notdir $<)

grid_collocate_miniapp.x: grid_collocate_miniapp.o $(ALL_OBJECTS)
	nvcc $(NVFLAGS) -o $@ $^ $(LIBS)

grid_collocate_unittest.x: grid_collocate_unittest.o $(ALL_OBJECTS)
	nvcc $(NVFLAGS) -o $@ $^ $(LIBS)

#EOF
