c234567
      subroutine init_cond(pos_src, pos_targ)

	real * 8 pos_src(3, 4)
	real * 8 pos_targ(3, 3)

	real * 8 xmax
	
	xmax = 0.5

	pos_src(1,1) = -xmax 
	pos_src(2,1) = -xmax
	
	pos_src(1,2) = xmax
	pos_src(2,2) = -xmax	

	pos_src(1,3) = -xmax
	pos_src(2,3) = xmax

	pos_src(1,4) = xmax
	pos_src(2,4) = xmax

	pos_targ(3,1) = -1.0
	pos_targ(3,3) = 1.0
      end subroutine init_cond

c234567 =============================================================
c       This subroutine takes in a complex vector (3 x N) and set it 
c       to 0 in parallel

      subroutine set_vector(vector, N, val, nthread)

c      include 'omp_lib.h'
      use omp_lib

      implicit none

      complex *16 vector(3, N)

      integer N, nthread

      complex * 16 val, uu

      integer i, j, tid, max_threads
   
      call omp_set_num_threads(nthread)


C$OMP PARALLEL SHARED (N, val, vector) PRIVATE(tid, i,j)
      tid = omp_get_thread_num()
      max_threads = omp_get_num_threads()
      
C$OMP DO SCHEDULE(DYNAMIC)
      do i=1, N
	do j=1,3
	  vector(j,i) = val
	enddo
      enddo
C$OMP END DO
      
c$OMP END PARALLEL

      return
      end subroutine set_vector
c234567
      program main

	implicit none

	integer nsource, ntarget, nthread, uu
	parameter(nsource = 4, ntarget = 3)

	complex * 16 fldsrc(3,nsource), fldsrc_near(3,nsource)
	complex * 16 fldsrc_far(3, nsource)

	complex * 16 fldtarg(3, ntarget), fldtarg_near(3, ntarget)
	complex * 16 fldtarg_far(3, ntarget)

	real * 8 pos_src(3, nsource), pos_targ(3, ntarget)

	complex * 16 val

	! set initial condition
	call init_cond(pos_src, pos_targ)

	val = 1.0
	nthread = 2

	call set_vector(fldsrc, nsource, val, nthread)
	call zero_vector(fldsrc, nsource, nthread)
	uu=1
      end program main

      