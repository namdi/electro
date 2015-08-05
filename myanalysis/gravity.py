import numpy as np

# general data directory
data_dir = "/home/namdi/Documents/School/UNC/Parallel_Time/Data"

path1 = data_dir + "/Results/2014_02/11"
ref_dir = path1 + "/sdirect_s256"
fname = ref_dir + "/t_0/charge.txt"

# load charge data
charge_src = np.loadtxt(fname).view('c16')

# make all charges positive
charge_src.real = abs(charge_src.real)

# save charges
fname = "/home/namdi/Desktop/charge.txt"
np.savetxt( fname, np.column_stack( [charge_src.real, charge_src.imag] ) )