# The purpose of this code is to have instructions to send src code from local computer to Killdevil

MY_BASE_DIR="/home/namdi/Documents/School/UNC/Parallel_Time"

# get the PyPFASST directory
FPFASST_DIR=${MY_BASE_DIR}/Code/fpfasst

# get the lustre scratch-space directory
lustre_dir=/lustre/scr/b/r/brandonn

# get the net scratch-space directory
nas02_dir=/nas02/home/b/r/brandonn

# get the killdevil scratch-space directory
kd_dir=$nas02_dir/fpfasst

dir=/electro

#To transfer a directory FROM my computer TO killdevil (Be on the local computer):
scp -r ${FPFASST_DIR}${dir} brandonn@killdevil.unc.edu:${kd_dir}