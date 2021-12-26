#include <stdio.h>
#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
    int dims[2];
    int periodic[2];
    int reorder = 1, q = 4, ndims = 2, maxdims = 2;
    int coordinates[2];
    int my_grid_rank;
    int coords[2];
    
    MPI_Comm grid_comm;
    MPI_Comm row_comm;

    dims[0] = dims[1] = q;
    periodic[0] = periodic[1] = 1;
    coords[0] = 0; coords[1] = 1;

    MPI_Init(&argc, &argv);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, reorder, &grid_comm);

    MPI_Comm_rank(grid_comm, &my_grid_rank);
    MPI_Cart_coords(grid_comm, my_grid_rank, maxdims, coordinates);
    //MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int* nneighbors);
    printf("Process rank %i has coordinates %i %i\n", my_grid_rank, coordinates[0], coordinates[1]);
    MPI_Finalize();
    return 0;
    

    MPI_Init(&argc, &argv);

    int index[] = { 4,1,1,1,1 };
    int edges[] = { 1,2,3,4,0,0,0,0 };
    MPI_Comm StarComm;
    MPI_Graph_create(MPI_COMM_WORLD, 5, index, edges, 1, &StarComm);

    int nneighbors;
    MPI_Graph_neighbors_count(StarComm, 1, &nneighbors);

    //int MPI_Graph_neighbors(MPI_Comm comm, int rank, int mneighbors, int* neighbors);

    cout << nneighbors << endl;

    MPI_Finalize();
    return 0;
}