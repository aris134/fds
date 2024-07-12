export OMP_TARGET_OFFLOAD=MANDATORY
mpifort -g -O3 -fiopenmp -fopenmp-targets=spir64 main_gpu_mpi_v7.f90 -o loop3d_mpi_gpu_v7