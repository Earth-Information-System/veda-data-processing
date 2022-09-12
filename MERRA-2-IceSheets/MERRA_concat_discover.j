#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o container_out_%04A_%04a
#SBATCH -e container_err_%04A_%04a

# This should be a one and done script to run the MERRA_concat_discover.py
# code. Run using:
# >>> sbatch MERRA_concat_discover.j
# (I did not test the log functionality yet.)

module load python/GEOSpyD/Min4.11.0_py3.9

srun -N1 -n1 -c1 python MERRA_concat_discover.py 


