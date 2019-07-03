#!/bin/bash
#SBATCH --job-name=pyatoa_test
#SBATCH --time=00:05:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263    
#SBATCH --clusters=maui_ancil         
#SBATCH --partition=nesi_prepost
#SBATCH --export=None
#SBATCH --output=test.out

# Kinda hacky but parameters are taken from positional arguments
event_id="$1"
model_num="$2"
step_count="$3"
cur_dir="$4"
work_dir="$5"
out_dir="$6"
suffix="$7"

# load anaconda and hdf5, start conda environment
echo "loading modules"
module load Anaconda3/5.2.0-GCC-7.1.0
module load HDF5/1.10.1-GCC-7.1.0 
source activate tomo

python /scale_wlg_persistent/filesets/project/nesi00263/PyPackages/pyatoa/pyatoa/scripts/seisflows/process.py \
--event_id ${event_id} \
--model_number ${model_num} \
--step_count ${step_count} \
--current_dir ${cur_dir} \
--working_dir ${work_dir} \
--output_dir ${out_dir} \
--suffix ${suffix}

echo "finished `date`"
