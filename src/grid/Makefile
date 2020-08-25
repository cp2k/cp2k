CFLAGS   = -fopenmp -g -O3 -std=c99 -march=native -Wall -Wextra -Wpedantic

all: grid_collocate_miniapp.x grid_collocate_unittest.x

grid_collocate_miniapp.o: grid_collocate_replay.h common/grid_library.h
grid_collocate_unittest.o: grid_collocate_replay.h common/grid_library.h
grid_collocate_replay.o: grid_collocate_replay.h grid_collocate.h grid_task_list.h common/grid_common.h common/grid_basis_set.h
grid_task_list.o: grid_task_list.h common/grid_common.h common/grid_constants.h common/grid_library.h ref/grid_ref_task_list.h cpu/grid_cpu_task_list.h
grid_collocate.o: grid_collocate.h ref/grid_ref_collocate.h

ref/grid_ref_task_list.o: ref/grid_ref_task_list.h ref/grid_ref_collocate.h common/grid_common.h common/grid_basis_set.h
ref/grid_ref_collocate.o: ref/grid_ref_collocate.h ref/grid_ref_prepare_pab.h common/grid_common.h common/grid_library.h
ref/grid_ref_prepare_pab.o: ref/grid_ref_prepare_pab.h common/grid_common.h common/grid_constants.h

cpu/grid_cpu_task_list.o: cpu/grid_cpu_task_list.h common/grid_basis_set.h

common/grid_library.o: common/grid_library.h common/grid_constants.h
common/grid_basis_set.o: common/grid_basis_set.h

%.o: %.c
	cd $(dir $<); $(CC) -c $(CFLAGS) $(notdir $<)

grid_collocate_miniapp.x: grid_collocate_miniapp.o \
                          grid_collocate_replay.o grid_task_list.o grid_collocate.o \
                          common/grid_library.o common/grid_basis_set.o \
                          cpu/grid_cpu_task_list.o \
                          ref/grid_ref_task_list.o ref/grid_ref_collocate.o ref/grid_ref_prepare_pab.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

grid_collocate_unittest.x: grid_collocate_unittest.o \
                           grid_collocate_replay.o grid_task_list.o grid_collocate.o \
                           common/grid_library.o common/grid_basis_set.o \
                           cpu/grid_cpu_task_list.o \
                           ref/grid_ref_task_list.o ref/grid_ref_collocate.o ref/grid_ref_prepare_pab.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -fv grid_*.o */grid_*.o grid_*.x

#EOF
