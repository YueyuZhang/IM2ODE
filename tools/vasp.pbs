#!/bin/bash -x

#PBS -l nodes=3:ppn=4,cput=08:10:00
#PBS -j oe
#PBS -q q1.3
#PBS -N zyy

#define MPI PATH
OMPI_HOME=/home/software/mpi/openmpi1.4.2-intel

# Setup the OpenMPI topology
n_proc=$(cat $PBS_NODEFILE | wc -l)

contexts=`~/bin/get_psm_sharedcontexts_max.sh`
 if [ "$?" = "0" ] ; then
  export PSM_SHAREDCONTEXTS_MAX=$contexts
 fi

cd $PBS_O_WORKDIR

cp INCAR_1 INCAR
python writekp.py 0.08
$OMPI_HOME/bin/mpirun --mca btl openib,self -machinefile $PBS_NODEFILE -np $n_proc ${HOME}/bin/v5211complex.omp &>log.out

cp CONTCAR POSCAR
cp INCAR_2 INCAR
python writekp.py 0.08
$OMPI_HOME/bin/mpirun --mca btl openib,self -machinefile $PBS_NODEFILE -np $n_proc ${HOME}/bin/v5211complex.omp &>log.out

cp CONTCAR POSCAR
cp INCAR_3 INCAR
python writekp.py 0.03
$OMPI_HOME/bin/mpirun --mca btl openib,self -machinefile $PBS_NODEFILE -np $n_proc ${HOME}/zhangyueyu/bin/vasp.5.optic &>log.out

exit 0
