#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <iomanip>
#include <stdlib.h>

#define MAX_VALUE INT_MAX

using namespace std;

double MinInArr(int n, int* Arr) {
    int min_n = MAX_VALUE;
    for (int i = 0; i < n; ++i) {
        min_n = min(min_n, Arr[i]);
    }
    return min_n;
}

double MinInArrParallel(int size, int n, int* Arr) {
    int n_proc_rank = n / size;
    int* array_proc_rank = new int[n_proc_rank];
    MPI_Scatter(Arr, n_proc_rank, MPI_INT, array_proc_rank, n_proc_rank, MPI_INT, 0, MPI_COMM_WORLD);

    /*for (int i = 0; i < n_proc_rank; ++i) {
        cout << array_proc_rank[i] << " ";
    }
    cout << endl;
    cout << endl;*/

    int min_value = MAX_VALUE;

    for (int i = 0; i < n_proc_rank; ++i) {
        min_value = min(min_value, array_proc_rank[i]);
    }

    int total_min = MAX_VALUE;
    MPI_Reduce(&min_value, &total_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    delete[] array_proc_rank;
    return total_min;
}

int main(int argc, char* argv[]) {
    int size, rank;
    double start, end, start_parallel, end_parallel;
    double time, time_parallel;
    int total_min, total_min_parallel;

    int* Arr = nullptr;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        cout << endl;
        cout << endl;
        cout << "  number of threads = " << size << endl;
        cout << "  -------------------------------------------------------------" << endl;
        cout << "                           Min In Arr                          " << endl;
        cout << "  -------------------------------------------------------------" << endl;
        cout << "         n  |       times        |    s   |      results       " << endl;
        cout << "  -------------------------------------------------------------" << endl;
    }
    for (int n = 100000; n <= 409600000; n *= 2) {
        if (rank == 0) {
            Arr = new int[n];
            for (int i = 0; i < n; ++i) {
                Arr[i] = rand();
            }

            start = MPI_Wtime();
            total_min = MinInArr(n, Arr);
            end = MPI_Wtime();
            time = end - start;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        start_parallel = MPI_Wtime();
        total_min_parallel = MinInArrParallel(size, n, Arr);
        MPI_Barrier(MPI_COMM_WORLD);
        end_parallel = MPI_Wtime();
        time_parallel = end_parallel - start_parallel;

        if (rank == 0) {
            cout << setw(10) << right << n << "  |  " \
                << fixed << setprecision(4) << time << "    " \
                << fixed << setprecision(4) << time_parallel << "  |  " \
                << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
                << fixed << setprecision(0) << total_min << "    " \
                << fixed << setprecision(0) << total_min_parallel << endl;
        }
    }
    delete[] Arr;

    MPI_Finalize();
    return 0;
}