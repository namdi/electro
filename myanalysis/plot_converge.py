#
# Plot error-vs-time for KS equation to demonstrate that PFASST is
# converging properly.
#
# This assumes that you have already run:
#
# $ ./main.exe probin.ks.ref.nml
# $ mpiexec -n 5 ./main.exe probin.ks.pf.nml
#
# Make sure that the LIBPFASST directory is in your PYTHONPATH:
#
# export PYTHONPATH=/path/to/libpfasst
#

import pylab
import pf

path = "/home/namdi/Documents/School/UNC/Parallel_Time/Data"

reference   = pf.io.read_avail(path + "/ref")
approximate = pf.io.read_avail(path + "/test")
errs, steps, iters, levels = pf.convergence.errors(reference, approximate)

pf.convergence.plot(errs, steps, iters, levels, color='k')

pylab.show()
