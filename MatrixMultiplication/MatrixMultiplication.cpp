#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>

#define n 150

using namespace std;

/*print arr*/
void ArrPrint(double(&Arr)[n][n]) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << Arr[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int main(int argc, char* argv[])
{
    int	size, rank;
    double start_time, end_time, start_time1, end_time1;
    double a[n][n] = { 0 };
    double b[n][n] = { 0 };
    double c[n][n] = { 0 };
    double c1[n][n] = { 0 };
    int portion, low_bound, upper_bound;

    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
    {
        srand(time(0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                a[i][j] = i + j;
                b[i][j] = i * j;
            }
        }

        //ArrPrint(a);
        //ArrPrint(b);

        start_time = MPI_Wtime();

        for (int i = 1; i < size; ++i)
        {
            portion = n / size - 1;
            low_bound = (i - 1) * portion;
            if ((i + 1 == size) && (n % (size - 1) != 0)) {
                upper_bound = n;
            }
            else {
                upper_bound = low_bound + portion;
            }
            MPI_Isend(&low_bound, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &request);
            MPI_Isend(&upper_bound, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &request);
            MPI_Isend(&a[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &request);
        }
    }

    MPI_Bcast(b, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank > 0)
    {
        MPI_Recv(&low_bound, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&a[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);

        for (int i = low_bound; i < upper_bound; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        MPI_Isend(&low_bound, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &request);
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &request);
        MPI_Isend(&c[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &request);
    }

    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Recv(&low_bound, 1, MPI_INT, i, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_bound, 1, MPI_INT, i, 5, MPI_COMM_WORLD, &status);
            MPI_Recv(&c[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
        }
        end_time = MPI_Wtime();
        cout << end_time - start_time << endl;

        //ArrPrint(c);

        start_time1 = MPI_Wtime();

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                c1[i][j] = 0;
                for (int k = 0; k < n; ++k) {
                    c1[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        end_time1 = MPI_Wtime();
        cout << end_time1 - start_time1 << endl;

        //ArrPrint(c1);

    }

    MPI_Finalize();
    return 0;
}