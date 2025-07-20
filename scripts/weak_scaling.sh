#!/bin/bash

echo "Submitting Weak Scaling jobs..."

TASKS_PER_NODE=8
CPUS_PER_TASK=14
LOCAL_SIZE=2048 # single node grid size
N_STEPS=200


for nodes in 1 2 4 8 16; do
    TOTAL_TASKS=$((nodes * TASKS_PER_NODE))
    
    case ${TOTAL_TASKS} in
        8) PX=4; PY=2; ;;
        16) PX=4; PY=4; ;;
        32) PX=8; PY=4; ;;
        64) PX=8; PY=8; ;;
        128) PX=16; PY=8; ;;
        *) echo "Invalid number of tasks: ${TOTAL_TASKS}"; exit 1; ;;
    esac

    GRID_SIZE_X=$((PX * LOCAL_SIZE))
    GRID_SIZE_Y=$((PY * LOCAL_SIZE))
    PARAMS="-x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -n ${N_STEPS}"

    JOB_NAME="weak_scale_${nodes}n_${TOTAL_TASKS}t"
    echo "Submitting job: ${JOB_NAME} with ${nodes} nodes and ${TOTAL_TASKS} tasks..."

    sbatch --nodes=${nodes} \
           --ntasks-per-node=${TASKS_PER_NODE} \
           --cpus-per-task=${CPUS_PER_TASK} \
           --job-name=${JOB_NAME} \
           --export=ALL,PROGRAM_ARGS="${PARAMS}", TEST_TYPE="weak" \
           go_dcgp.sbatch
done

echo "All Weak Scaling jobs submitted."