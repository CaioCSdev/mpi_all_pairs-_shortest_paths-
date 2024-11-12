#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define tag 0

typedef struct {
  int n;
  int *vals;
} Matrix;

typedef struct {
  int procs;         // number of processes
  MPI_Comm comm;     // communicator for the grid
  MPI_Comm row_comm; // row communicator
  MPI_Comm col_comm; // column communicator
  int q;             // dimension of the grid
  int my_row;        // my row number
  int my_col;        // my column number
  int my_rank;       // my rank in the grid
} GRID;

Matrix read_file(char *filename);
void print_location();
void print_matrix(Matrix matrix);
int min_in_ij(Matrix *a, Matrix *b, int i, int j);
void calculate(Matrix *A, Matrix *B, Matrix *C);
int usable_procs(int num_procs, int n);

void Fox(int n, GRID *grid, Matrix *local_A, Matrix *local_B, Matrix *local_C) {
  Matrix *temp_A;
  int step;
  int bcast_root;
  int n_bar;
  int source;
  int dest;
  MPI_Status status;

  n_bar = n / grid->q;
  for (int i = 0; i < n_bar * n_bar; i++)
    local_C->vals[i] = INT_MAX;

  /* Calculate addresses for circular shift of B */
  source = (grid->my_row + 1) % grid->q;
  dest = (grid->my_row + grid->q - 1) % grid->q;

  temp_A = (Matrix *)malloc(sizeof(Matrix));
  temp_A->n = n_bar;
  temp_A->vals = (int *)malloc(n_bar * n_bar * sizeof(int));

  // printf("row: %d, col: %d\n", grid->my_row, grid->my_col);
  for (step = 0; step < n; step++) {
    bcast_root = (grid->my_row + step) % grid->q;
    if (bcast_root == grid->my_col) {
      MPI_Bcast(local_A->vals, local_A->n, MPI_INT, bcast_root, grid->row_comm);
      calculate(local_A, local_B, local_C);
    } else {
      MPI_Bcast(temp_A->vals, temp_A->n, MPI_INT, bcast_root, grid->row_comm);
      calculate(temp_A, local_B, local_C);
    }
    MPI_Send(local_B->vals, local_B->n, MPI_INT, dest, tag, grid->col_comm);
    MPI_Recv(local_B->vals, local_B->n, MPI_INT, source, tag, grid->col_comm,
             &status);
  }
  // print_matrix(*local_C);
}

int main(int argc, char **argv) {
  int n, my_rank, num_procs;
  int dimensions[2], coordinates[2];
  int periods[2] = {1, 1};
  int varying_coords[2] = {0, 1};
  Matrix matrix;
  GRID grid;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  if (my_rank == 0) {
    matrix = read_file(argv[1]);
    n = matrix.n;
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (my_rank != 0) {
    matrix.n = n;
    matrix.vals = (int *)malloc(n * n * sizeof(int));
  }

  MPI_Bcast(matrix.vals, n * n, MPI_INT, 0, MPI_COMM_WORLD);

  grid.procs = usable_procs(num_procs, matrix.n);
  grid.q = sqrt(grid.procs);
  dimensions[0] = dimensions[1] = grid.q;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &grid.comm);
  if (grid.comm == MPI_COMM_NULL) {
    /*Cannot use this processor, the grid cannot comport it*/
    MPI_Finalize();
    return 0;
  }
  MPI_Comm_rank(grid.comm, &grid.my_rank);
  MPI_Cart_coords(grid.comm, grid.my_rank, 2, coordinates);
  grid.my_row = coordinates[0];
  grid.my_col = coordinates[1];
  MPI_Cart_sub(grid.comm, varying_coords, &grid.row_comm);
  varying_coords[0] = 1;
  varying_coords[1] = 0;
  MPI_Cart_sub(grid.comm, varying_coords, &grid.col_comm);
  /*
  The code above creates a grid of 2 dimensions divided in q x q processes
  and assigns a row and column communicator to each process.
  */

  Matrix local_A, local_B, local_C;
  local_A.n = local_B.n = local_C.n = matrix.n / grid.q;
  local_A.vals = (int *)malloc(local_A.n * local_A.n * sizeof(int));
  local_B.vals = (int *)malloc(local_B.n * local_B.n * sizeof(int));
  local_C.vals = (int *)malloc(local_C.n * local_C.n * sizeof(int));

  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n; j++) {
      int row = i / local_A.n;
      int col = j / local_A.n;
      if (grid.my_row == row && grid.my_col == col) {
        int local_i = i % local_A.n;
        int local_j = j % local_A.n;
        local_A.vals[local_j + local_i * local_A.n] =
            matrix.vals[j + i * matrix.n];
        local_B.vals[local_j + local_i * local_B.n] =
            matrix.vals[j + i * matrix.n];
      }
    }
  }

  for (int i = 0; i < matrix.n - 2; i++) {

    Fox(matrix.n, &grid, &local_A, &local_B, &local_C);
  local_A = local_C;
  local_B = local_C;
  }

  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n; j++) {
      int row = i / local_A.n;
      int col = j / local_A.n;
      if (grid.my_row == row && grid.my_col == col) {
        int local_i = i % local_A.n;
        int local_j = j % local_A.n;

        matrix.vals[j + i * matrix.n] =
            local_C.vals[local_j + local_i * local_C.n];
      }
    }
  }

  if (my_rank == 0) {
    print_matrix(matrix);
  }
  // if (my_rank == 1) {
  //   printf("my_rank is: %d, local_A n is: %d, procs: %d\n", my_rank,
  //   local_A.n,
  //          num_procs);
  //   print_matrix(local_A);
  // }

  // if (my_rank == 0) {
  // for (int i = 0; i < matrix.n - 2; i++)
  //   calculate(&matrix, &matrix, &matrix);
  //   print_matrix(matrix);
  // }
  // print_location();

  MPI_Finalize();
  return 0;
}

int usable_procs(int num_procs, int n) {
  /* returns the largest p which p <= num_procs, q*q = p and n % q == 0 */
  for (int q = sqrt(num_procs); q > 0; q--) {
    if (n % q == 0)
      return q * q;
  }

  return 1;
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

  matrix.vals = (int *)malloc(matrix.n * matrix.n * sizeof(int));
  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n; j++) {
      fscanf(file, "%d", &matrix.vals[j + i * matrix.n]);
      if (i != j && matrix.vals[j + i * matrix.n] == 0)
        matrix.vals[j + i * matrix.n] = INT_MAX;
    }
  }

  fclose(file);

  return matrix;
}

int min_in_ij(Matrix *a, Matrix *b, int i, int j) {
  int min = a->vals[j + i * a->n];
  int candidate;

  for (int k = 0; k < a->n; k++) {
    if (a->vals[k + i * a->n] == INT_MAX || b->vals[j + k * a->n] == INT_MAX)
      continue;
    candidate = a->vals[k + i * a->n] + b->vals[j + k * a->n];
    if (candidate < min)
      min = candidate;
  }

  return min;
}

void calculate(Matrix *a, Matrix *b, Matrix *c) {
  for (int i = 0; i < a->n; i++) {
    for (int j = 0; j < a->n; j++) {
      int old_min = c->vals[j + i * a->n];
      int new_min = min_in_ij(a, b, i, j);
      if (new_min < old_min)
        c->vals[j + i * a->n] = new_min;
    }
  }
}

void print_matrix(Matrix matrix) {
  int val;

  for (int i = 0; i < matrix.n; i++) {
    for (int j = 0; j < matrix.n; j++) {
      if (matrix.vals[j + i * matrix.n] == INT_MAX)
        val = 0;
      else
        val = matrix.vals[j + i * matrix.n];

      if (j == matrix.n - 1) {
        printf("%d", val);
      } else {
        printf("%d ", val);
      }
    }
    printf("\n");
  }
}

void print_location() {
  pid_t pid = getpid();
  char hostname[HOST_NAME_MAX + 1];

  gethostname(hostname, HOST_NAME_MAX + 1);
  printf("Hi from Hostname: %s PID: %d\n", hostname, pid);
}
