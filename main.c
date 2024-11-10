#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
  printf("some values: #%d -> %s \n", argc, argv[0]);
  MPI_Init(&argc, &argv);



  MPI_Finalize();
  return 0;
}
