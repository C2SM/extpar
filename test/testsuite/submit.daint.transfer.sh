#!/bin/bash -l
#
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --partition=xfer
#SBATCH --output="transfer.log"

command="cp"

origin=/store/c2sm/extpar_raw_data/linked_data
destination=input-data

mkdir -p $destination

for file in $(cat transfer.txt); do
    echo "copy $origin/$file to $destination"
    srun -n $SLURM_NTASKS $command $origin/$file $destination/.
done
