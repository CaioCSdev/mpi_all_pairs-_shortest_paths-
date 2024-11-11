#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct {
  int n;
  int **vals;
} Matrix;

Matrix read_file(char *filename);
void print_location();
void print_matrix(Matrix matrix);

int main(int argc, char **argv) {
  int my_rank;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    Matrix matrix = read_file(argv[1]);
    print_matrix(matrix);
  }

  // print_location();

  MPI_Finalize();
  return 0;
}

Matrix read_file(char *filename) {
  FILE *file = fopen(filename, "r");

  if (file == NULL) {
    printf("Error: Could not open file %s\n", filename);
    Matrix empty_matrix = {0, NULL};
    return empty_matrix;
  }

  Matrix matrix;
  fscanf(file, "%d", &matrix.n);

  matrix.vals = (int **)malloc(matrix.n * sizeof(int *));
  for (int i = 0; i < matrix.n; i++) {
    matrix.vals[i] = (int *)malloc(matrix.n * sizeof(int));
    for (int j = 0; j < matrix.n; j++) {
      fscanf(file, "%d", &matrix.vals[i][j]);
    }
  }

  fclose(file);

  return matrix;
}

void print_matrix(Matrix matrix) {
  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n; j++)
      printf("%d ", matrix.vals[i][j]);
    printf("\n");
  }
}

void print_location() {
  pid_t pid = getpid();
  char hostname[HOST_NAME_MAX + 1];

  gethostname(hostname, HOST_NAME_MAX + 1);
  printf("Hi from Hostname: %s PID: %d\n", hostname, pid);
}
