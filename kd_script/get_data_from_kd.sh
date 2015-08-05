# The purpose of this code is to get data FROM killdevil TO my local computer

DATA_DIR=/home/namdi/Documents/School/UNC/Parallel_Time/Data

# get the lustre scratch-space directory
lustre_dir=/lustre/scr/b/r/brandonn

# get the netscr scratch-space
netscr_dir=/netscr/brandonn

# kill devil data directory
kd_data_dir=${lustre_dir}/Data

# the source of the dta to be copied
dat_source=/sdirect_s256_f9

# the kill devil directory of the data to be copied
kd_dir=${kd_data_dir}${dat_source}

# the local directory in which to store the info
local_dir=${DATA_DIR}/kd

# TO make a transfer of a directory FROM Killdevil to my local computer
scp -r brandonn@killdevil.unc.edu:${kd_dir} ${local_dir}
