#!/bin/bash

echo "Submitting OpenMP scaling tests"

NODES=1
TASKS_PER_NODE=1
GRID_SIZE="-x 16384 -y 16384" #2^14 x 2^14
N_STEPS="-n 200"

#14 and 28 corresponds to 2 numa regions
for threads in 1 2 4 8 14 28 56 84 112; do
    JOB_NAME="omp_scale_${threads}t"
    echo "Submitting job: ${JOB_NAME} with ${threads} threads..."

    sbatch --nodes=${NODES} \
           --ntasks-per-node=${TASKS_PER_NODE} \
           --cpus-per-task=${threads} \
           --job-name=${JOB_NAME} \
           --export=ALL,PROGRAM_ARGS="${GRID_SIZE} ${N_STEPS}", TEST_TYPE="omp" \
           go_dcgp.sbatch
done

echo "All OpenMp jobs submitted"