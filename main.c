#include <limits.h>
#include "main.h"
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  print_location();

  MPI_Finalize();
  return 0;
}

void print_location()
{
  pid_t pid = getpid();
  char hostname[HOST_NAME_MAX + 1];

  gethostname(hostname, HOST_NAME_MAX + 1);
  printf("Hi from Hostname: %s PID: %d\n", hostname, pid);
}
