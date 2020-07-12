#$ -cwd -V
#$ -l nodes=1,ppn=16,tpp=1
#$ -l h_rt=48:00:00
#$ -m be
#$ -M eepsh@leeds.ac.uk
#$ -N modellingHTI

#!/bin/bash
export OMP_NUM_THREADS=16
./fdmodeling
