# Hello_Parallel_MPI
:white_check_mark:

# Pi_Monte_Carlo_method
:white_check_mark:

ssh USER@ant.apmath.spbu.ru

job.sh =============================

#!/bin/sh

module load openmpi

mpirun ~/Mpi/Hello/hello


#!/bin/sh

module load openmpi

for((i = 2; i <= 16; i++))

do

echo $i

mpirun -np $((i)) ~/Mpi/Pi_Monte_Carlo/Pi_Monte_Carlo_method

done


=====================================

chmod +x job.sh

sbatch -p gnu job.sh

Makefile =============================

hello: Hello.cpp
        mpicc Hello.cpp -o hello
clean:
        rm -f hello
