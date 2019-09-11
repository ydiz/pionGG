#!/bin/bash

export LD_LIBRARY_PATH=/home/yidizhao/cooley/src/boost_1_70_0/stage/lib:/home/yidizhao/cooley/install/gcc-6.4/lib:/home/yidizhao/cooley/install/gcc-6.4/lib64:/home/yidizhao/cooley/install/hdf5/lib:$LD_LIBRARY_PATH

NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 1))

mpirun -f $COBALT_NODEFILE -n $PROCS ./jackknife --mpi 1.1.1.1 --ensemble @ensemble@ --target @target@ --file_p3 @file_p3@ --time_cutoff_end @time_cutoff_end@ > @outfile@

