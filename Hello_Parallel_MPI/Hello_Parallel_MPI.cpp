#include <stdio.h>
#include <mpi.h>

void PrintHello(int rank, int numtasks) {
    printf("Hello MPI from process = %2d, total number of processes: %d\n", rank, numtasks);
}

int main(int argc, char** argv) {

    int numtasks, rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);


    PrintHello(4, 4);
    

    MPI_Finalize();
    return 0;
}