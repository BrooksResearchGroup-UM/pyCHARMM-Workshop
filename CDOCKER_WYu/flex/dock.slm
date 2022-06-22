#!/bin/bash

#Created by Yujin Wu (wyujin@umich.edu)
#at 2020/04/10

#SBATCH -A brooks
#SBATCH -p gpu2080 
#SBATCH --job-name=flex
#SBATCH --time=5:00:00
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --output=slurm.out

############################################################################################
## Prepare environment 
############################################################################################

## Module load in satyr
module load pycharmm/0.5

## Docking
python standard.py

 
