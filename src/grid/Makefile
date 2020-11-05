.PHONY : all clean

all: grid_collocate_miniapp.x grid_collocate_unittest.x

clean:
	rm -fv grid_*.o */grid_*.o grid_*.x

CFLAGS := -fopenmp -g -O3 -std=c99 -march=native -Wall -Wextra -Wpedantic

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
               ref/grid_ref_prepare_pab.o

%.o: %.c $(ALL_HEADERS)
	cd $(dir $<); $(CC) -c $(CFLAGS) $(notdir $<)

grid_collocate_miniapp.x: grid_collocate_miniapp.o $(ALL_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lblas

grid_collocate_unittest.x: grid_collocate_unittest.o $(ALL_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lblas

#EOF
