CFLAGS   = -g -O3 -std=c99 -march=native -Wall -Wextra

all: grid_collocate_miniapp.x

grid_collocate_cpu.o: grid_prepare_pab.o
grid_collocate_replay.o: grid_collocate_cpu.o

%.o: %.c
	$(CC) -c $(CFLAGS) $<

grid_collocate_miniapp.x: grid_collocate_miniapp.o grid_collocate_replay.o grid_collocate_cpu.o grid_prepare_pab.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -fv grid_*.o grid_*.x

#EOF
