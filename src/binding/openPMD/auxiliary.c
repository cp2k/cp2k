#include "auxiliary.h"

void CP2K_MPI_Comm_f2c(MPI_Fint comm_f, MPI_Comm *comm_c) {
  *comm_c = MPI_Comm_f2c(comm_f);
}