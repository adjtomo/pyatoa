EVENT_ID=2015p768477
BASE_DIR='/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/tests/seisflows_test/hikurangi_trial/twoevent'
python ./process.py \
--event_id ${EVENT_ID} \
--model_number m00 \
--step_count 0 \
--current_dir ${BASE_DIR}/scratch/solver/${EVENT_ID} \
--output_dir ${BASE_DIR}/pyatoa.output \
--working_dir ${BASE_DIR} \
--suffix test


