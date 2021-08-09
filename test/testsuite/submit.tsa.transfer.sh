#!/bin/bash -l
#
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --partition=debug
#SBATCH --output="transfer.log"

command="cp"

origin=/store/c2sm/extpar_raw_data/linked_data
destination=input-data

mkdir -p $destination

for file in $(cat transfer.txt); do
    echo "copy $origin/$file to $destination"
    $command $origin/$file $destination/.
done
