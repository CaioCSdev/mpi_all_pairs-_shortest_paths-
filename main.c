#include <limits.h>
#include <math.h>
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

int min_in_ij(Matrix matrix, int i, int j) {
  int min = matrix.vals[i][j];
  int candidate;

  for (int k = 0; k < matrix.n; k++) {
    if (matrix.vals[i][k] == INT_MAX || matrix.vals[k][j] == INT_MAX)
      continue;
    candidate = matrix.vals[i][k] + matrix.vals[k][j];
    if (candidate < min)
      min = candidate;
  }

  return min;
}

void step(Matrix matrix) {
  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n; j++) {
      matrix.vals[i][j] = min_in_ij(matrix, i, j);
    }
  }
}

int main(int argc, char **argv) {
  int my_rank;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    Matrix matrix = read_file(argv[1]);
    for (int i = 0; i < matrix.n - 2; i++)
      step(matrix);
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
      if (i != j && matrix.vals[i][j] == 0)
        matrix.vals[i][j] = INT_MAX;
    }
  }

  fclose(file);

  return matrix;
}

void print_matrix(Matrix matrix) {
  int val;

  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n - 1; j++) {
      val = matrix.vals[i][j] == INT_MAX ? 0 : matrix.vals[i][j];
      printf("%d ", val);
    }
    val = matrix.vals[i][matrix.n - 1] == INT_MAX
      ? 0
      : matrix.vals[i][matrix.n - 1];
    printf("%d", val);
    printf("\n");
  }
}

void print_location() {
  pid_t pid = getpid();
  char hostname[HOST_NAME_MAX + 1];

  gethostname(hostname, HOST_NAME_MAX + 1);
  printf("Hi from Hostname: %s PID: %d\n", hostname, pid);
}
