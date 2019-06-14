# positional arguments
# event id, model number, step_count, current_dir, working_dir, output_dir, suffix
EVENT_ID=2015p768477
MODEL='m00'
STEP_COUNT=0
BASE_DIR='/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/tests/seisflows_test/hikurangi_trial/twoevent'
sbatch ./run.sh \
${EVENT_ID} \
${MODEL} \
${STEP_COUNT} \
${BASE_DIR}/scratch/solver/${EVENT_ID} \
${BASE_DIR} \
${BASE_DIR}/pyatoa.output \
test


