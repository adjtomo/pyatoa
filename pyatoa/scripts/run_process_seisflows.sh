#!/bin/bash
#SBATCH --job-name=pyatoa_test
#SBATCH --time=00:10:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263    
#SBATCH --clusters=maui_ancil         
#SBATCH --partition=nesi_prepost
#SBATCH --export=None

# Kinda hacky but parameters are taken from positional arguments
event_id="$1"
model_num="$2"
syn_dir="$3"
work_dir="$4"
out_dir="$5"

# load anaconda and hdf5, start conda environment
echo "loading modules"
module load Anaconda3/5.2.0-GCC-7.1.0
module load HDF5/1.10.1-GCC-7.1.0 
source activate tomo

python /scale_wlg_persistent/filesets/project/nesi00263/PyPackages/pyatoa/pyatoa/scripts/process_seisflows.py \
-i ${event_id} -m ${model_num} -p ${syn_dir} -w ${work_dir} -o ${out_dir}


echo "finished `date`"
