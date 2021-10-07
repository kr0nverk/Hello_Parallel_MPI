#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>

using namespace std;

int main(int argc, char* argv[])
{
    int	size, rank, rows, offset, offset2, source;

    double start_time_cons, end_time_cons, start_time, end_time;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const unsigned n = 81;

    /*
    double** a = new double* [n];
    double** b = new double* [n];
    double** c = new double* [n];
    for (int i = 0; i < n; ++i) {
        a[i] = new double[n];
        b[i] = new double[n];
        c[i] = new double[n];
    }
    */
    double a[n][n], b[n][n], c[n][n];
    int num_workers = size - 1;

    if (rank == 0)
    {

        cout << "Init" << endl;
        srand((unsigned)(time(0)));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                a[i][j] = (double)rand();
                b[i][j] = (double)rand();
                //cout << a[i][j] << " ";
            }
            //cout << endl;
        }

        start_time_cons = MPI_Wtime();

        //rows = n / num_workers;
        offset = 0;

        cout << "1" << endl;

        for (int dest = 1; dest <= num_workers; ++dest)
        {
            rows = n / num_workers;
            offset = (dest - 1) * rows;
            if ((dest + 1 == size) && (n % num_workers != 0)) {
                offset2 = n;
            }
            else {
                offset2 = offset + rows;
            }
            MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&a[offset][0], rows * n, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&b, n * n, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            offset = offset + rows;
        }

        cout << "2" << endl;
        for (int i = 1; i <= num_workers; ++i)
        {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&c[offset][0], rows * n, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
        }
        cout << "3" << endl;

        end_time_cons = MPI_Wtime();
        cout << end_time_cons - start_time_cons << endl;
    }


    if (rank > 0)
    {
        source = 0;
        MPI_Recv(&offset, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&a, rows * n, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&b, n * n, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);

        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                c[i][k] = 0.0;
                for (int j = 0; j < n; ++j) {
                    c[i][k] = c[i][k] + a[i][j] * b[j][k];
                }
            }
        }

        MPI_Send(&offset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&c, rows * n, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
    /*
    for (int i = 0; i < n; ++i) {
        delete[] a[i];
        delete[] b[i];
        delete[] c[i];
    }
    delete[]a;
    delete[]b;
    delete[]c;

    */


    MPI_Finalize();
    return 0;
}