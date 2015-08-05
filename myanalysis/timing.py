import numpy as np

"""
This data is from 2014_05/solver_timings

the fmm tolerance based on iprec flag.
    if( iprec .eq. -2 ) epsfmm=.5d-0 
    if( iprec .eq. -1 ) epsfmm=.5d-1
    if( iprec .eq. 0 ) epsfmm=.5d-2
    if( iprec .eq. 1 ) epsfmm=.5d-3
    if( iprec .eq. 2 ) epsfmm=.5d-6
    if( iprec .eq. 3 ) epsfmm=.5d-9
    if( iprec .eq. 4 ) epsfmm=.5d-12
    if( iprec .eq. 5 ) epsfmm=.5d-15
    if( iprec .eq. 6 ) epsfmm=0
"""

# number of processor
Nproc = 64

# maximum iterations
niter = 6

# serial multirate (7 iterations, omp 2)
smr3 = 3821.

# serial FMM (7 iterations, omp 2)
sfmm3 = 4820.

"""
# parallel FMM 7 iterations
pfmm03 = 274.
pfmm13 = 313.
pfmm23 = 510.
pfmm33 = 710.

# parallel MRFMM 7 iterations
pmr03 = 156.
pmr13 = 209.
pmr23 = 427.
pmr33 = 522.
"""

# parallel FMM 6 iterations
pfmmn23 = 167.
pfmmn13 = 299.
pfmm03 = 200.
pfmm13 = 260.
pfmm23 = 492.
pfmm33 = 687.

# parallel MRFMM 6 iterations
pmrn23 = 124.
pmrn13 = 127.
pmr03 = 148.
pmr13 = 192.
pmr23 = 414.
pmr33 = 510.

# array of parallel timings
pfmm = np.array((pfmm03, pfmm13, pfmm23, pfmm33))

pmr = np.array((pmr03, pmr13, pmr23, pmr33))

# speedup 
# serial solutions converged in 6 iterations
speedup_fmm = sfmm3*(6./7.)/pfmm
speedup_mr = sfmm3*(6./7.)/pmr

# efficiency vs the 
eff_fmm = speedup_fmm/Nproc
eff_mr = speedup_mr/Nproc
eff_mr_fmm = eff_mr / eff_fmm

print "speedup: fmm, mr"
print speedup_fmm
print speedup_mr
print "\nefficiency: fmm, mr"
print eff_fmm
print eff_mr
print "\nefficiency/speedup: fmm/mr"
print eff_mr_fmm
print
print 
print "runtimes: FMM, MR"
print pfmm
print pmr
