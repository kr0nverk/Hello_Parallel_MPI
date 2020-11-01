#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
	int proc_rank, proc_num;
	double x_rand, y_rand, start_time_cons, end_time_cons, start_time, end_time;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		start_time_cons = MPI_Wtime();

		int n = 16000000;
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

		end_time_cons = MPI_Wtime();

		printf("Points: %d\n", n);
		printf("Pi: %10.5f\n", monte_carlo_pi);
		printf("Time: %5.3f\n", end_time_cons - start_time_cons);

	}


	MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime();
	
	int n = 16000000;
	int rank_hit = 0;
	srand((unsigned)(time(0)));

	int n_proc_rank = n / proc_num;
	for (int i = 0; i < n_proc_rank; ++i) {
		x_rand = (double)rand() / (double)RAND_MAX;
		y_rand = (double)rand() / (double)RAND_MAX;
		if (x_rand * x_rand + y_rand * y_rand <= 1.0) {
			rank_hit++;
		}
	}
	int hit = 0;
	MPI_Reduce(&rank_hit, &hit, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	double monte_carlo_pi = hit * 4 / (double)n;

	MPI_Barrier(MPI_COMM_WORLD);
	end_time = MPI_Wtime();

	if (proc_rank == 0) {
		printf("Points: %d\n",n);
		printf("Pi: %10.5f\n", monte_carlo_pi);
		printf("Time: %5.3f\n", end_time - start_time);
	}
	
	MPI_Finalize();
	return 0;
}