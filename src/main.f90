!
! Simple example of using LIBPFASST.
! This is based off of the main function from mpi-ndarray
!
! To run the pfasst code:
!
! read the integer command line parameters in the following order:
! nthreads, solver, accuracy, iterations, nlevs, fnodes, indir, outdir 

! To run:
! mpirun -n 8./main.exe 4 $FMM 3 5 2 5 $MY_DATA_DIR/test/t_0 $MY_DATA_DIR/test
!
!	MPI threads: 8	OpenMP threads/ MPI node: 4	solver: $FMM
! 	accuracy: 3	iterations: 5	nlevs: 2	fine nodes: 5
!	input directory: $MY_DATA_DIR/test/t_0		output directory: $MY_DATA_DIR/test
!
!===============================================
! Function declarations
!===============================================

!
! parse command line
!

subroutine parse_cmd_line(nthreads, accuracy, solver, iterations, nlevs, fnodes, indir, outdir)

  implicit none

  integer, intent(out) :: nthreads, accuracy, solver, iterations, nlevs, fnodes

  character(len=*), intent(out) :: indir

  character(len=*), intent(out) :: outdir

  character(len=200) :: arg

  ! read the integer command line parameters in the following order:
  ! nthreads, solver, accuracy, iterations, nlevs, fnodes, indir, outdir 

  ! read nthreads
  call get_command_argument(1, arg)
  read(arg, *) nthreads

  ! read the solver
  call get_command_argument(2, arg)
  read(arg, *) solver

  ! read the accuracy
  call get_command_argument(3, arg)
  read(arg,*) accuracy

  ! read the number of iterations
  call get_command_argument(4, arg)
  read(arg, *) iterations

  call get_command_argument(5, arg)
  read(arg, *) nlevs

  call get_command_argument(6, arg)
  read(arg, *) fnodes

  ! read the indir
  call get_command_argument(7, indir)

  ! read the outdir if it is defined
  if(iargc() .ge. 6) then
    call get_command_argument(8, outdir)
  endif
  
 ! the value of accuracy corresponds to the following error tolerances for the FMM potential.
 ! For the E-field subtract 3 digits of precision.
 ! -2 => tolerance =.5d0
 ! -1 => tolerance =.5d-1
 ! 0 => tolerance =.5d-2
 ! 1 => tolerance =.5d-3
 ! 2 => tolerance =.5d-6
 ! 3 => tolerance =.5d-9
 ! 4 => tolerance =.5d-12
 ! 5 => tolerance =.5d-15
  
end subroutine parse_cmd_line

!===============================================
! Run
!===============================================

program main
  use pfasst
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use pf_mod_ndarray
  !use verlet
  !use encap ! part of verlet
  use feval
  use hooks
  use transfer
  use my_utils
  !use probin

  implicit none

  !----------------------------------------------
  ! Interface
  !----------------------------------------------
  interface
    subroutine parse_cmd_line(nthreads, accuracy, solver, iterations, nlevs, fnodes, indir, outdir)

      integer, intent(out) :: nthreads, accuracy, solver, iterations, nlevs, fnodes

      character(len=*), intent(out) :: indir

      character(len=*), intent(out) :: outdir

      character(len=200) :: arg
    end subroutine parse_cmd_line
  end interface

  !integer, parameter 		:: nlevs = 2

  type(pf_pfasst_t) 		:: pf
  type(pf_comm_t)   		:: comm

  type(ndarray),      target 	:: q1
  type(pf_sweeper_t), target 	:: sweeper
  type(pf_encap_t),   target 	:: encapsulation

  integer			:: err, l, nlevs, fnodes, nsources, ntargets, & ! originally: nvars(levs), nnodes(nlevs)
				   nsteps, accuracy, solver, nthreads, iterations, N, length, &
				   distribution, dim, ifcharge, ifdipole, ifsource, iftarget, &
				   ifpot, iffld, i

  integer, allocatable		:: nvars(:), nnodes(:)

  double precision		:: dt
  character(len=300)		:: fname, fname2, indir, outdir

  complex*16, pointer 		:: charge_src(:), dipstr_src(:)

  real*8, pointer		:: dipvec_src(:,:)

  character(len=100)		:: solver_char

  !
  ! initialize mpi
  !

  call mpi_init(err)
  if (err .ne. 0) &
       stop "ERROR: Can't initialize MPI."

  !
  ! read command line / load parameters
  !

  outdir = ""

  call parse_cmd_line(nthreads, accuracy, solver, iterations, nlevs, fnodes, indir, outdir)
  !call probin_init('/home/namdi/Documents/School/UNC/Parallel_Time/Code/fpfasst/PFASST/examples/electro/probin.nml')  

  ! for the debugger
  !nthreads 	= 4
  !solver 	= 3
  !accuracy 	= 3
  !iterations 	= 4
  !nlevs		= 2
  !fnodes	= 5

  !indir = "/home/namdi/Documents/School/UNC/Parallel_Time/Data/test/t_0"
  !outdir = "/home/namdi/Documents/School/UNC/Parallel_Time/Data/test"
  
  allocate(nvars(nlevs))
  allocate(nnodes(nlevs))

  ! store the solver in a character variable
  select case(solver)
    case(SOLVE_DIRECT)
      solver_char = "SOLVE_DIRECT"
    case(SOLVE_FMM)
      solver_char = "SOLVE_FMM"
    case(SOLVE_FMM_NEAR)
      solver_char = "SOLVE_FMM_NEAR"
    case(SOLVE_FMM_NEAR_FAR)
      solver_char = "SOLVE_FMM_NEAR_FAR"
    case default
      stop "STOP: invalid solver chosen"
  end select

  length = len_trim(indir)
  fname = indir(1:length) // "/info.txt"

  ! load parameters
  call load_parameters(fname, nsources, ntargets, distribution, dt, nsteps)

  N = nsources + ntargets

  ! load charge and dipoles

  allocate( charge_src(nsources) )
  allocate( dipstr_src(nsources) )
  allocate( dipvec_src(3, nsources) )

  fname = indir(1:length) // "/charge.txt"
  call mu_load_charge(fname, charge_src, nsources)

  fname = indir(1:length) // "/dipstr.txt"
  fname2 = indir(1:length) // "/dipvec.txt"
  call mu_load_dipole(fname, fname2, dipstr_src, dipvec_src, nsources)

  ! Solver parameters
  
  ifpot		= 0
  iffld		= 1

  ifcharge	= 1
  ifdipole	= 0

  ifsource	= 1
  iftarget	= 0 

  ! 
  ! pfasst parameters
  !
  
  ! set the number of variables per level
  do i=1, nlevs
    nvars(i) = N   
  enddo

  ! set the number of nodes per level (coarsest to finest)
  if (nlevs == 1) then
    nnodes(nlevs) = fnodes
  else if (nlevs > 1) then    
    do i=1, nlevs    
      nnodes(i) = 2**i + 1
    enddo
  endif  

  dim		= 3	! dimension in shape

  !
  ! Write output file about pfasst parameters
  !

  !
  ! initialize pfasst using 2 levels
  !

  call ndarray_encap_create(encapsulation)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_explicit_create(sweeper, eval_f1)
  !call verlet_create(sweeper, acceleration, hamiltonian, 0.0_pfdp, 0.0_pfdp)
  call pf_pfasst_create(pf, comm, nlevs)

  pf%niters = iterations
  pf%qtype  = SDC_GAUSS_LOBATTO ! + SDC_PROPER_NODES
  pf%window = PF_WINDOW_BLOCK
 
  !pf%echo_timings = .true.

  ! Namdi: I added this data member  
  pf%ctx_char = outdir  

  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     allocate(pf%levels(l)%shape(3))

     pf%levels(l)%shape  = [ 3, nvars(l), 2]

     pf%levels(l)%nvars  = product(pf%levels(l)%shape)
     pf%levels(l)%nnodes = nnodes(l)

     call feval_create_workspace( pf%levels(l)%ctx,  nsources, ntargets, nthreads, accuracy, &
				  ifcharge, ifdipole, ifsource, iftarget, solver, ifpot, iffld, l, &
				  charge_src, dipstr_src, dipvec_src, nlevs )

     pf%levels(l)%encap       => encapsulation
     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%sweeper     => sweeper
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  if (pf%rank == 0) then
     print *, 'nlevs:  ', nlevs
     print *, 'nvars:  ', pf%levels(:)%nvars
     print *, 'nnodes: ', pf%levels(:)%nnodes
     print *, 'output: ', trim(pf%ctx_char)
     length = len_trim(pf%ctx_char)
     print *, 'len_trim(outdir): ', length
     print *, ""
     print *, "nthreads: ", nthreads
     print *, "solver: ", trim(solver_char)
     print *, "accuracy: ", accuracy
     print *, "iterations: ", iterations
  end if

  !
  ! run
  !

  call ndarray_create_simple(q1, pf%levels(nlevs)%shape)
  call initial(indir, q1, nsources, ntargets)

  if (len_trim(outdir) > 0) then  
     call dump_mkdir(outdir, len_trim(outdir))
     call pf_add_hook(pf, nlevs, PF_POST_SWEEP, my_dump_hook)
     !call pf_add_hook(pf, nlevs, PF_POST_SWEEP, dump_hook)     
  end if

  call mpi_barrier(MPI_COMM_WORLD, err)

  call pf_pfasst_run(pf, c_loc(q1), dt, tend=0.d0, nsteps=nsteps)  
  
  !
  ! cleanup
  !

  deallocate(q1%flatarray)
  
 
  do l = 1, nlevs
     call feval_destroy_workspace(pf%levels(l)%ctx)
  end do

  deallocate(nvars)
  deallocate(nnodes)
  deallocate(charge_src)
  deallocate(dipvec_src)
  deallocate(dipstr_src)

  ! verlet
  !call verlet_destroy(sweeper)

  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(err)  

end program main
