#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]) {

    int numtasks, rank;
    int n = 10000;
    double* blocklens = new double[n];
    double* indices = new double[n];

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_numtasks);

    if (proc_rank == 0) {

        for (int i = 0, i < n; ++i) {
            blocklens[i] = n - i;
            indices[i] = i * n + i;
        }
        MPI_Type_indexed(n, blocklens, indices, &UTMatrixType, &ElemType);

        
        MPI_Send(&)

    }
    else {
        
    }




    delete[] blocklens;
    delete[] indices;


    MPI_Finalize();
}