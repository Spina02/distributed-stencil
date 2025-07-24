#!/bin/sh  
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH -t 00:10:00
#SBATCH --job-name=serial_stencil
##SBATCH --exclusive
#SBATCH --mem=0

EXEC=build/serial

# =======================================================
module purge
# module load gcc/12.2.0

${EXEC} -x 16384 -y 16384 -n 500 -s