# The purpose of this code is to have instructions to send the Data directory from local computer to Killdevil

# local data directory
DATA_DIR=/home/namdi/Documents/School/UNC/Parallel_Time/Data

# the directory with the data we want to load
fdir=/sfmm_s256_a4_f9
#fdir=/t_0

dir_name=${DATA_DIR}/kd${fdir}

kd_dir=/lustre/scr/b/r/brandonn/Data${fdir}

#To transfer a directory FROM my computer TO killdevil (Be on the local computer):
scp -r ${dir_name} brandonn@killdevil.unc.edu:${kd_dir}
