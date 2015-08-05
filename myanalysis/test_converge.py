"""
Plot error-vs-time for Laplace equation equation to demonstrate that PFASST is
converging properly.

This assumes that you have already have made a reference solution and a non-reference solution

In KS example, make sure you run:
$ ./main.exe probin.ks.ref.nml
$ mpiexec -n 5 ./main.exe probin.ks.pf.nml

Make sure that the LIBPFASST directory is in your PYTHONPATH:

export PYTHONPATH=/path/to/libpfasst
"""

import pylab
import pf

path1 = "/home/namdi/Documents/School/UNC/Parallel_Time/Data"
path2 = "/home/namdi/Documents/School/UNC/Parallel_Time/Data/Results/2013_11/01"
path3 = "/home/namdi/Documents/School/UNC/Parallel_Time/Data/Results/2013_11/04"
path4 = "/home/namdi/Documents/School/UNC/Parallel_Time/Data/kd"

reference   = pf.io.read_avail(path4 + "/test")
approximate = pf.io.read_avail(path4 + "/test")

#errs, steps, iters, levels = pf.convergence.errors(reference, approximate)

pf.my_test.errors(reference, approximate)

#pf.convergence.plot(errs, steps, iters, levels, color='k')

#pylab.show()
