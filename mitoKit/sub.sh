#!/usr/bin/sh
#SBATCH -p BP10
#SBATCH -J MITO-KIT
#SBATCH -N 1
#SBATCH -n 4
./mito-kit.pl
