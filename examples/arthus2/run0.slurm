#!/bin/sh
#SBATCH -J cluster
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-core=1
#SBATCH --mail-user=Ronan.Guivarch@enseeiht.fr
#SBATCH --mail-type=ALL

module purge
module load intel/14.0.2.144 intelmpi/4.1.3.049 
module load mpiprof_intelmpi/intelmpi4_compilo14
module list

cd $SLURM_SUBMIT_DIR 

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun /users/p0989/guivarch/Cluster/bin/clusters param0.in > out.${SLURM_JOBID}
ls cluster.final.* | wc -l > nbclusters

