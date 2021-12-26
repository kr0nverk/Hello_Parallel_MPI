#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <iomanip>
#include <stdlib.h>

using namespace std;

/*MonteCarloPi*/
double MonteCarloPi(int n) {
	double x_rand, y_rand;
	int rank_hit = 0;
	srand((unsigned)(time(0)));

	for (int i = 0; i < n; ++i) {
		x_rand = (double)rand() / (double)RAND_MAX;
		y_rand = (double)rand() / (double)RAND_MAX;
		if (x_rand * x_rand + y_rand * y_rand <= 1.0) {
			rank_hit++;
		}
	}
	double monte_carlo_pi = rank_hit * 4 / (double)n;

	return monte_carlo_pi;
}

/*MonteCarloPiParallel*/
double MonteCarloPiParallel(int n, int size) {
	double x_rand, y_rand;
	int rank_hit = 0;
	srand((unsigned)(time(0)));

	int n_proc_rank = n / size;
	for (int i = 0; i < n_proc_rank; ++i) {
		x_rand = (double)rand() / (double)RAND_MAX;
		y_rand = (double)rand() / (double)RAND_MAX;
		if (x_rand * x_rand + y_rand * y_rand <= 1.0) {
			rank_hit++;
		}
	}
	int hit = 0;
	MPI_Reduce(&rank_hit, &hit, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	double monte_carlo_pi_parallel = hit * 4 / (double)n;

	return monte_carlo_pi_parallel;
}
int main(int argc, char* argv[]) {
	int rank, size;
	double monte_carlo_pi, monte_carlo_pi_parallel, start, end, start_parallel, end_parallel;
	double time, time_parallel;

	int n = 16000000;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (rank == 0) {
		cout << endl;
		cout << endl;
		cout << "  number of threads = " << size << endl;
		cout << "  -------------------------------------------------------------" << endl;
		cout << "                         Pi Monte Carlo                        " << endl;
		cout << "  -------------------------------------------------------------" << endl;
		cout << "         n  |       times        |    s   |      results       " << endl;
		cout << "  -------------------------------------------------------------" << endl;
	}

	for (int n = 1600000; n < 160000000; n += 16000000) {

		if (rank == 0) {
			start = MPI_Wtime();
			monte_carlo_pi = MonteCarloPi(n);
			end = MPI_Wtime();
			time = end - start;
			/*
			cout << "Points: " << n << endl;
			cout << "Pi: " << monte_carlo_pi << endl;
			cout << "Time: " << time << endl;
			*/
		}

		MPI_Barrier(MPI_COMM_WORLD);
		start_parallel = MPI_Wtime();
		monte_carlo_pi_parallel = MonteCarloPiParallel(n, size);
		MPI_Barrier(MPI_COMM_WORLD);
		end_parallel = MPI_Wtime();
		time_parallel = end_parallel - start_parallel;

		if (rank == 0) {
			/*
			cout << "Points: " << n << endl;
			cout << "Pi: " << monte_carlo_pi_parallel << endl;
			cout << "Time: " << time_parallel << endl;
			*/

			cout << setw(10) << right << n << "  |  " \
				<< fixed << setprecision(4) << time << "    " \
				<< fixed << setprecision(4) << time_parallel << "  |  " \
				<< fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
				<< fixed << setprecision(4) << monte_carlo_pi << "    " \
				<< fixed << setprecision(4) << monte_carlo_pi_parallel << endl;

		}
	}
	MPI_Finalize();
	return 0;
}