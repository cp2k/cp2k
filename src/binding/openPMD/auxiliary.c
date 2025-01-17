#include "auxiliary.h"

#include <string.h>

void CP2K_MPI_Comm_f2c(MPI_Fint comm_f, MPI_Comm *comm_c) {
  *comm_c = MPI_Comm_f2c(comm_f);
}

size_t CP2K_strlen(char const *s) { return strlen(s); }