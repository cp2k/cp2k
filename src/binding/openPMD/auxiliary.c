#include "auxiliary.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void CP2K_MPI_Comm_f2c(MPI_Fint comm_f, MPI_Comm *comm_c) {
  *comm_c = MPI_Comm_f2c(comm_f);
}

size_t CP2K_strlen(char const *s) { return strlen(s); }

void CP2K_free(void *ptr) { free(ptr); }

void CP2K_log(char const *s) {
  FILE *file = fopen("log.txt", "a");

  if (file == NULL) {
    perror("Error opening log.txt for appending");
    return;
  }
  fprintf(file, "%s\n", s);
  fclose(file);
}
