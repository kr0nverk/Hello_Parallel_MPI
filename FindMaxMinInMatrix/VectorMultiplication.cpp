#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <iomanip>


using namespace std;

double VectorMultiplication(int n, int* A, int* B) {
    long long vector_multiplication = 0;

    for (int i = 0; i < n; ++i) {
        vector_multiplication += (long long)A[i] * B[i];
    }

    return vector_multiplication;
}

double VectorMultiplicationParallel(int size, int n, int* A, int* B) {
    long long vector_multiplication_parallel = 0;

    int n_proc_rank = n / size;
    int* A_n_proc = new int[n_proc_rank];
    int* B_n_proc = new int[n_proc_rank];

    MPI_Scatter(A, n_proc_rank, MPI_INT, A_n_proc, n_proc_rank, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, n_proc_rank, MPI_INT, B_n_proc, n_proc_rank, MPI_INT, 0, MPI_COMM_WORLD);

    long long local_multi = 0;
    for (int i = 0; i < n_proc_rank; ++i) {
        local_multi += (long long)A_n_proc[i] * B_n_proc[i];
    }

    MPI_Reduce(&local_multi, &vector_multiplication_parallel, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    delete[] A_n_proc;
    delete[] B_n_proc;

    return vector_multiplication_parallel;
}

int main(int argc, char* argv[]) {
    int size, rank;
    double start, end, start_parallel, end_parallel;
    double time, time_parallel;
    long long vector_multiplication, vector_multiplication_parallel;
    
    int* A = nullptr;
    int* B = nullptr;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        cout << endl;
        cout << endl;
        cout << "  number of threads = " << size << endl;
        cout << "  -------------------------------------------------------------" << endl;
        cout << "                     Vector Multiplication                     " << endl;
        cout << "  -------------------------------------------------------------" << endl;
        cout << "         n  |       times        |    s   |      results       " << endl;
        cout << "  -------------------------------------------------------------" << endl;
    }

    for (int n = 625000; n <= 160000000; n *= 2) {
        if (rank == 0) {
            A = new int[n];
            B = new int[n];
            for (int i = 0; i < n; ++i) {
                A[i] = rand() % 1000;
                B[i] = rand() % 1000;
                //cout << A[i] << " " << B[i] << endl;
            }
            //cout << endl;

            start = MPI_Wtime();
            vector_multiplication = VectorMultiplication(n, A, B);
            end = MPI_Wtime();
            time = end - start;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        start_parallel = MPI_Wtime();
        vector_multiplication_parallel = VectorMultiplicationParallel(size, n, A, B);
        MPI_Barrier(MPI_COMM_WORLD);
        end_parallel = MPI_Wtime();
        time_parallel = end_parallel - start_parallel;


        if (rank == 0) {
            cout << setw(10) << right << n << "  |  " \
                << fixed << setprecision(4) << time << "    " \
                << fixed << setprecision(4) << time_parallel << "  |  " \
                << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
                << fixed << setprecision(0) << vector_multiplication << "    " \
                << fixed << setprecision(0) << vector_multiplication_parallel << endl;
        }
    }
    delete[] A;
    delete[] B;

    MPI_Finalize();
    return 0;
}


/*#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>

#define n 6

using namespace std;

/*print arr
void ArrPrint(double Arr[n * n]) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << Arr[i*n+j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int main(int argc, char* argv[])
{
    int	size, rank;
    double start_time, end_time, start_time1, end_time1;
    /*double a[n][n] = { 0 };
    double b[n][n] = { 0 };
    double c[n][n] = { 0 };
    double c1[n][n] = { 0 };

    //int n = 4;
    double* a = new double[n * n];
    double* b = new double[n * n];
    double* c = new double[n * n];
    double* c1 = new double[n * n];

    /* int *ary = new int[sizeX*sizeY]
    ary[i][j] is then rewritten as
        ary[i * sizeY + j]

    for (int i = 0; i < n; ++i) {
        for (int j = 0; i < n; ++i) {
            /*a[i][j] = 0;
            b[i][j] = 0;
            c[i][j] = 0;
            c1[i][j] = 0;
            a[i * n + j] = 0;
            b[i * n + j] = 0;
            c[i * n + j] = 0;
            c1[i * n + j] = 0;
        }
    }

    int portion, portion_ost, low_bound, upper_bound;


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
                /*a[i][j] = i+j;
                b[i][j] = i*j;
                a[i * n + j] = i + j;
                b[i * n + j] = i * j;
            }
        }

        ArrPrint(a);
        ArrPrint(b);

        start_time = MPI_Wtime();

        for (int i = 1; i < size; ++i)
        {
            /*
            portion = n / size - 1;
            low_bound = (i - 1) * portion;
            if ((i + 1 == size) && (n % (size - 1) != 0)) {
                upper_bound = n;
            }
            else {
                upper_bound = low_bound + portion;
            }
            portion = (n / (size - 1)) * n;
            low_bound = (i - 1) * portion;
            if ( ( n % (size - 1) != 0 ) && ( i + 1 == size ) ) {
                upper_bound = n;
            }
            else {
                upper_bound = low_bound + portion;
            }

            //cout << portion << " " << low_bound << " " << endl;
            //MPI_Isend(&low_bound, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &request);
            //MPI_Isend(&upper_bound, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &request);
            //MPI_Isend(&a[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &request);
            MPI_Send(&low_bound, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&upper_bound, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
            //MPI_Send(&a[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
            MPI_Send(&(a[low_bound]), (upper_bound - low_bound) * n, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
        }
    }

    MPI_Bcast(b, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank > 0)
    {
        MPI_Recv(&low_bound, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
        //MPI_Recv(&a[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
        MPI_Recv(&(a[low_bound]), (upper_bound - low_bound) * n, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
        
        for (int i = low_bound; i < upper_bound; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    /*c[i][j] += a[i][k] * b[k][j];
                    c[i * n + j] += a[i * n + k] * b[k * n + j];
                }
            }
        }

        //MPI_Isend(&low_bound, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &request);
        //MPI_Isend(&upper_bound, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &request);
        //MPI_Isend(&c[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &request);
        MPI_Send(&low_bound, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
        MPI_Send(&upper_bound, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);
        //MPI_Send(&c[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
        MPI_Send(&(c[low_bound]), (upper_bound - low_bound) * n, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Recv(&low_bound, 1, MPI_INT, i, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_bound, 1, MPI_INT, i, 5, MPI_COMM_WORLD, &status);
            //MPI_Recv(&c[low_bound][0], (upper_bound - low_bound) * n, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
            MPI_Recv(&(c[low_bound]), (upper_bound - low_bound)* n, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
        }
        end_time = MPI_Wtime();
        cout << end_time - start_time << endl;
        
        ArrPrint(c);

        start_time1 = MPI_Wtime();

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                //c1[i][j] = 0;
                c1[i * n + j] = 0;
                for (int k = 0; k < n; ++k) {
                    //c1[i][j] += a[i][k] * b[k][j];
                    c1[i * n + j] += a[i * n + k] * b[k * n + j];
                }
                //if (c[i][j] != c1[i][j]) {
                if (c[i*n+j] != c1[i*n+j]) {
                    cout << "alert";
                }
            }
        }
        end_time1 = MPI_Wtime();
        cout << end_time1 - start_time1 << endl;

        ArrPrint(c1);

    }

    /*
    for (int i = 0; i < n; ++i) {
        delete[] a[i];
    }
    delete[]a;

    for (int i = 0; i < n; ++i) {
        delete[] b[i];
    }
    delete[]b;

    for (int i = 0; i < n; ++i) {
        delete[] c[i];
    }
    delete[]c;

    for (int i = 0; i < n; ++i) {
        delete[] c1[i];
    }
    delete[]c1;
    MPI_Finalize();
   
    return 0;
}*/