!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use pf_mod_ndarray

  !use encap ! part of verlet
  use my_utils
  use my_type  
  use electro
  use tree
  use print_mod

  !use fortfmm_utils

  implicit none  

  ! advection work type
  type :: ad_work_t
     type(c_ptr) :: ffft, ifft				! FFT, inverse FFT
     complex(pfdp), pointer :: wk1(:)                   ! work space
     complex(pfdp), pointer :: ddx1(:), lap1(:), ks1(:) ! operators
  end type ad_work_t

contains

  !
  ! This function creates the workspace on the level, as in the electro object and trees
  ! construct electro object/ tree
  !

  subroutine feval_create_workspace(ctx, nsources, ntargets, nthreads, iprec, ifcharge, ifdipole, ifsource, iftarget, &
				    solver, ifpot, iffld, lev, charge_src, dipstr_src, dipvec_src, nlevs)

    type(c_ptr), intent(out) :: ctx

    integer,     intent(in)  :: nsources, ntargets, nthreads, iprec, ifcharge, ifdipole, ifsource, iftarget, &
				solver, ifpot, iffld, lev, nlevs

    complex*16, intent(in)   :: dipstr_src(nsources),  charge_src(nsources)

    real*8, intent(in)	     :: dipvec_src(3, nsources)

    type(electro_t), pointer :: work

    ! allocate work pointer
    allocate(work)
    ctx = c_loc(work)
    
    ! electro constructor
    call electro_construct(work, nsources, ntargets, nthreads, iprec, ifcharge, ifdipole, &
				ifsource, iftarget, solver, ifpot, iffld, lev, charge_src, dipstr_src, dipvec_src, nlevs)

  end subroutine feval_create_workspace

  !
  ! This function deallocates the memory
  !

  subroutine feval_destroy_workspace(ctx)

    type(c_ptr), intent(in) :: ctx

    type(electro_t), pointer :: work    

    call c_f_pointer(ctx, work)

    ! electro destructor
    call electro_destruct(work)

    deallocate(work)
    nullify(work)

  end subroutine feval_destroy_workspace

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initial(fdir, q0, nsources, ntargets)

    type(ndarray), intent(inout) :: q0

    integer, intent(in) :: nsources, ntargets

    character(len=*), intent(in) :: fdir

    integer :: length, N, error, i, j, k, idx

    double precision, allocatable :: q(:,:,:)
  
    N = nsources + ntargets

    ! allocate memory
    allocate( q(3, N, 2), stat=error )

    if(error .ne. 0) then
      print *, "error: couldn't allocate array q"
    endif
 
    ! load the initial condition in a multiple dimension manner
    call mu_load_ic(fdir, q, nsources, ntargets)

    ! Enter data into the flat array

    idx = 1
    do k=1,2
      do j=1, N
	do i=1,3
	  q0%flatarray(idx) = q(i,j,k)
	  idx = idx+1
	enddo
      enddo
    enddo
  
    deallocate(q, stat=error)
  
    if(error .ne. 0) then
      print *, "error: couldn't allocate q! "
    endif

  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, yex)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    integer    :: nvars, i, ii, nbox
    real(pfdp) :: tol, x

  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.

  subroutine eval_f1(yptr, t, level, ctx, f1ptr) !orginally had m, the mth node
    type(c_ptr), intent(in), value :: yptr, f1ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level
    !integer, intent(in), optional  :: m

    type(electro_t), pointer :: work
    real(pfdp),      pointer :: y(:,:,:), f1(:,:,:), pos_src(:,:), pos_targ(:,:)

    complex(pfdp), pointer	:: pot_src(:), pot_far_src(:), pot_near_src(:), &
				    pot_targ(:), pot_far_targ(:), pot_near_targ(:), &
				    field_src(:,:), field_near_src(:,:), field_far_src(:,:), &
				    field_targ(:,:),field_far_targ(:,:), field_near_targ(:,:), &
				    charge_src(:), charge_targ(:), dipstr_src(:)

    real(pfdp), pointer		:: dipvec_src(:,:), accel_src(:,:), accel_targ(:,:)

    integer 	:: nsources, ntargets, N, ifcharge, ifdipole, ifsource, iftarget, nthreads ,&
		    ifpot, ifpotsrc, ifpotsrc_far, ifpotsrc_near, &
		    ifpottarg, ifpottarg_far, ifpottarg_near, iffld, iffldsrc, iffldsrc_far, &
		    iffldsrc_near, iffldtarg, iffldtarg_far, iffldtarg_near, solver

    type(tree_t), pointer :: tree

    integer :: ier, node, i, j

    ! testing. don't need mpole's
    integer :: ifpotsrc_mpole, ifpottarg_mpole, iffldsrc_mpole, iffldtarg_mpole, ifnear_only
    complex*16, allocatable :: pot_mpole_src(:), pot_mpole_targ(:), field_mpole_src(:,:), field_mpole_targ(:,:)

    ! set the current time node
    ! Note: could have problems with things being called at m+1, before m see pf_explicit.f90

    !if ( present(m) ) then
      !node = m
    !else
      !node = -1
    !endif

    call c_f_pointer(ctx, work)

    ! 3D array [3 N, 2]
    y	=> array3(yptr)	
    f1 	=> array3(f1ptr)
    
    nthreads	= work%nthreads

    nsources	= work%nsources
    ntargets	= work%ntargets
    N		= work%N
    solver	= work%solver

    ifsource	= work%ifsource
    iftarget	= work%iftarget

    ifcharge	= work%ifcharge
    ifdipole	= work%ifdipole

    ifpotsrc		= work%ifpotsrc
    ifpotsrc_far	= work%ifpotsrc_far
    ifpotsrc_near	= work%ifpotsrc_near

    ifpottarg		= work%ifpottarg
    ifpottarg_far	= work%ifpottarg_far
    ifpottarg_near	= work%ifpottarg_near

    iffldsrc		= work%iffldsrc
    iffldsrc_far	= work%iffldsrc_far
    iffldsrc_near	= work%iffldsrc_near

    iffldtarg		= work%iffldtarg
    iffldtarg_far	= work%iffldtarg_far
    iffldtarg_near	= work%iffldtarg_near
    
    charge_src	=> work%charge_src
    charge_targ => work%charge_targ

    dipstr_src	=> work%dipstr_src
    dipvec_src	=> work%dipvec_src

    pos_src	=> y(1:3, 1:nsources, 1)
    pos_targ	=> y(1:3, nsources+1:N, 1)    

    pot_src		=> work%pot_src
    pot_far_src		=> work%pot_far_src
    pot_near_src	=> work%pot_near_src

    pot_targ		=> work%pot_targ
    pot_far_targ	=> work%pot_far_targ
    pot_near_targ	=> work%pot_near_targ

    field_src		=> work%field_src
    field_far_src	=> work%field_far_src
    field_near_src	=> work%field_near_src

    field_targ		=> work%field_targ
    field_far_targ	=> work%field_far_targ
    field_near_targ	=> work%field_near_targ

    accel_src		=> work%accel_src
    accel_targ		=> work%accel_targ

    tree		=>work%tree

    ifnear_only 	= work%ifnear_only

    !--------------------------------------
    ! Testing
    !--------------------------------------
    
    allocate( pot_mpole_src(nsources) )
    allocate( pot_mpole_targ(ntargets) )
    allocate( field_mpole_src(3, nsources) )
    allocate( field_mpole_targ(3, ntargets) )
    
    ifpotsrc_mpole 	= 0
    ifpottarg_mpole	= 0
    iffldsrc_mpole	= 0
    iffldtarg_mpole	= 0    

    
    !------------------------------------------

    
    if (solver .eq. SOLVE_DIRECT) then
      
      call direct(nsources, pos_src, ifcharge,charge_src, ifdipole, dipstr_src, dipvec_src, & 
			   ifpotsrc, pot_src, iffldsrc, field_src, ntargets, pos_targ, ifpottarg, pot_targ, &
			  iffldtarg, field_targ, nthreads)           

      !call my_direct(field_src, pot_src, pos_src, nsources, charge_src, field_targ, pot_targ, &
	!		      pos_targ, ntargets, nthreads)
    
    else

      call tree_create(tree, pos_src, pos_targ)

      ! make the tree at the beginning of the timestep
      !if (node .eq. 1) then
	!call tree_create(tree, pos_src, pos_targ)
      !endif

      ! set parameters for FMM near-field only
      !if (solver .eq. SOLVE_FMM_NEAR) then

	!if (level .eq. 1) then ! coarse level

	  ! work%solver	= SOLVE_FMM_NEAR_FAR
	   !call set_solver(work)
	   !work%solver = SOLVE_FMM_NEAR

	!else if ( level .eq. 2. .and. (node .eq. 2 .or. node .eq. 4) ) then ! fine level

	  !call set_solver(work)
	  !ifnear_only = 1

	!endif 

      !endif !(solver .eq. SOLVE_FMM_NEAR)
    
      
      call fmm(ier, tree%iprec, nsources, pos_src, ifcharge, charge_src, ifdipole, dipstr_src, dipvec_src, &
		ifpotsrc, pot_src, iffldsrc, field_src, ntargets, pos_targ, ifpottarg, pot_targ, iffldtarg, &
		field_targ, tree%ladder, tree%nlev, tree%nboxes, tree%nbox, tree%epsfmm, tree%lused7, nthreads, tree%iisource, &
		tree%iwlists, tree%lwlists, tree%iitarget, tree%size, tree%wlists, tree%ntot, &
		ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
		iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
		ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
		iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
		pot_far_src,  field_far_src,  pot_mpole_src,  field_mpole_src,  pot_near_src,  field_near_src, &
		pot_far_targ, field_far_targ, pot_mpole_targ, field_mpole_targ, pot_near_targ, field_near_targ, &
		ifnear_only )

    endif ! (solver .eq. SOLVE_DIRECT)

    ! combine the near and far field for the respective FMM
    if ( solver .eq. SOLVE_FMM_NEAR .or. solver .eq. SOLVE_FMM_NEAR_FAR) then
      do i=1, nsources      
	!field_src(1, i) = field_near_src(1,i) + field_far_src(1,i)
	!field_src(2, i) = field_near_src(2,i) + field_far_src(2,i)
	!field_src(3, i) = field_near_src(3,i) + field_far_src(3,i)
	field_src(1:3,i) = field_near_src(1:3, i) + field_far_src(1:3, i)
      enddo

      do i=1, ntargets
	!field_targ(1, i) = field_near_targ(1,i) + field_far_targ(1,i)
	!field_targ(2, i) = field_near_targ(2,i) + field_far_targ(2,i)
	!field_targ(3, i) = field_near_targ(3,i) + field_far_targ(3,i)
	field_targ(1:3, i) = field_near_targ(1:3,i) + field_far_targ(1:3,i)
      enddo
      
    endif

    ! calculate the accelerations
    call calc_accel(accel_src, field_src, charge_src, nsources) 
    call calc_accel(accel_targ, field_targ, charge_targ, ntargets)

    ! store the velocity    
    f1(1:3, 1:N, 1) 		= y(1:3, 1:N, 2)

    ! store the acceleartion   
    f1(1:3, 1:nsources, 2) 	= accel_src(1:3, 1:nsources)
    f1(1:3, nsources+1:N, 2) 	= accel_targ(1:3, 1:ntargets)

    ! clean up: these values are used in testing

    deallocate( pot_mpole_src )
    deallocate( pot_mpole_targ )
    deallocate( field_mpole_src )
    deallocate( field_mpole_targ )
    !------------------------------------
  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_accel(accel, field, charge, N)

    integer, intent(in) 	:: N

    complex*16, intent(in)	:: field(3,N), charge(N)

    real*8, intent(out)		:: accel(3,N)

    integer ::  i,j 

    do j=1, N
      do i=1,3
	accel(i,j) = real( field(i,j) ) * real( charge(j) )
      enddo
    enddo 
    
  end subroutine

  !-----------------------------------------------
  ! Verlet functions
  !-----------------------------------------------

  subroutine hamiltonian(t, yptr, H)
    type(c_ptr), intent(in), value :: yptr
    real(pfdp),  intent(in)        :: t
    real(pfdp),  intent(out)       :: H
    
  end subroutine hamiltonian


!   subroutine acceleration(yptr, t, level, ctx, aptr)    
!     ! Namdi's variables
!     type(c_ptr), intent(in), value :: yptr, aptr, ctx
!     real(pfdp),  intent(in)        :: t
!     integer,     intent(in)        :: level
!     !integer, intent(in), optional  :: m


!     type(electro_t), pointer :: work
!     real(pfdp),      pointer :: y(:,:,:), f1(:,:,:), pos_src(:,:), pos_targ(:,:)

!     !Verlet
!     real(pfdp), pointer :: q(:, :), p(:, :), qdot(:, :), pdot(:, :) !Namdi: , m(:)

!     complex(pfdp), pointer	:: pot_src(:), pot_far_src(:), pot_near_src(:), &
! 				    pot_targ(:), pot_far_targ(:), pot_near_targ(:), &
! 				    field_src(:,:), field_near_src(:,:), field_far_src(:,:), &
! 				    field_targ(:,:),field_far_targ(:,:), field_near_targ(:,:), &
! 				    charge_src(:), charge_targ(:), dipstr_src(:)

!     real(pfdp), pointer		:: dipvec_src(:,:), accel_src(:,:), accel_targ(:,:)

!     integer 	:: nsources, ntargets, N, ifcharge, ifdipole, ifsource, iftarget, nthreads ,&
! 		    ifpot, ifpotsrc, ifpotsrc_far, ifpotsrc_near, &
! 		    ifpottarg, ifpottarg_far, ifpottarg_near, iffld, iffldsrc, iffldsrc_far, &
! 		    iffldsrc_near, iffldtarg, iffldtarg_far, iffldtarg_near, solver

!     type(tree_t), pointer :: tree

!     integer :: ier, node, i

!     ! testing. don't need mpole's
!     integer :: ifpotsrc_mpole, ifpottarg_mpole, iffldsrc_mpole, iffldtarg_mpole, ifnear_only
!     complex*16, allocatable :: pot_mpole_src(:), pot_mpole_targ(:), field_mpole_src(:,:), field_mpole_targ(:,:)   

!     !np = nparticles(yptr)
!     !nd = ndim(yptr)

!     !print *, "in accleration()"
!     !------------------------------------
!     ! VERLET
!     q => get_q(yptr); qdot => get_q(aptr)
!     p => get_p(yptr); pdot => get_p(aptr)
!     !m => get_m(yptr); 
!     !------------------------------------    
!     !print *, "after verlet pointer assign"
!     call c_f_pointer(ctx, work)

!     ! 3D array [3 N, 2]
!     !y	=> array3(yptr)	
!     !f1 	=> array3(aptr)
    
!     !print *, "after my verlet pointer assign"
!     nthreads	= work%nthreads

!     nsources	= work%nsources
!     ntargets	= work%ntargets
!     N		= work%N
!     solver	= work%solver

!     ifsource	= work%ifsource
!     iftarget	= work%iftarget

!     ifcharge	= work%ifcharge
!     ifdipole	= work%ifdipole

!     ifpotsrc		= work%ifpotsrc
!     ifpotsrc_far	= work%ifpotsrc_far
!     ifpotsrc_near	= work%ifpotsrc_near

!     ifpottarg		= work%ifpottarg
!     ifpottarg_far	= work%ifpottarg_far
!     ifpottarg_near	= work%ifpottarg_near

!     iffldsrc		= work%iffldsrc
!     iffldsrc_far	= work%iffldsrc_far
!     iffldsrc_near	= work%iffldsrc_near

!     iffldtarg		= work%iffldtarg
!     iffldtarg_far	= work%iffldtarg_far
!     iffldtarg_near	= work%iffldtarg_near
    
!     charge_src	=> work%charge_src
!     charge_targ => work%charge_targ

!     dipstr_src	=> work%dipstr_src
!     dipvec_src	=> work%dipvec_src

!     !pos_src	=> y(1:3, 1:nsources, 1)
!     !pos_targ	=> y(1:3, nsources+1:N, 1)    

!     pos_src	=> q(1:3, 1:nsources)
!     pos_targ	=> q(1:3, nsources+1:N)        
        
!     pot_src		=> work%pot_src
!     pot_far_src		=> work%pot_far_src
!     pot_near_src	=> work%pot_near_src

!     pot_targ		=> work%pot_targ
!     pot_far_targ	=> work%pot_far_targ
!     pot_near_targ	=> work%pot_near_targ

!     field_src		=> work%field_src
!     field_far_src	=> work%field_far_src
!     field_near_src	=> work%field_near_src

!     field_targ		=> work%field_targ
!     field_far_targ	=> work%field_far_targ
!     field_near_targ	=> work%field_near_targ

!     accel_src		=> work%accel_src
!     accel_targ		=> work%accel_targ

!     tree		=>work%tree

!     ifnear_only 	= work%ifnear_only

!     ! test: no problems
!     !--------------------------------------
!     ! Testing
!     !--------------------------------------
!     allocate( pot_mpole_src(nsources) )
!     allocate( pot_mpole_targ(ntargets) )
!     allocate( field_mpole_src(3, nsources) )
!     allocate( field_mpole_targ(3, ntargets) )
    
!     ifpotsrc_mpole 	= 0
!     ifpottarg_mpole	= 0
!     iffldsrc_mpole	= 0
!     iffldtarg_mpole	= 0    
    
!     !------------------------------------------

    
!     if (solver .eq. SOLVE_DIRECT) then

!       call direct(nsources, pos_src, ifcharge,charge_src, ifdipole, dipstr_src, dipvec_src, & 
! 			  ifpotsrc, pot_src, iffldsrc, field_src, ntargets, pos_targ, ifpottarg, pot_targ, &
! 			  iffldtarg, field_targ, nthreads)
!     else

!       ! make the tree at the beginning of the timestep
!       if (node .eq. 1) then
! 	call tree_create(tree, pos_src, pos_targ)
!       endif

!       ! set parameters for FMM near-field only
!       if (solver .eq. SOLVE_FMM_NEAR) then

! 	if (level .eq. 1) then ! coarse level

! 	   work%solver	= SOLVE_FMM_NEAR_FAR
! 	   call set_solver(work)
! 	   work%solver = SOLVE_FMM_NEAR

! 	else if ( level .eq. 2. .and. (node .eq. 2 .or. node .eq. 4) ) then ! fine level

! 	  call set_solver(work)
! 	  ifnear_only = 1

! 	endif 

!       endif !(solver .eq. SOLVE_FMM_NEAR)
    
      
!       call fmm(ier, tree%iprec, nsources, pos_src, ifcharge, charge_src, ifdipole, dipstr_src, dipvec_src, &
! 		ifpotsrc, pot_src, iffldsrc, field_src, ntargets, pos_targ, ifpottarg, pot_targ, iffldtarg, &
! 		field_targ, tree%ladder, tree%nlev, tree%nboxes, tree%nbox, tree%epsfmm, tree%lused7, nthreads, tree%iisource, &
! 		tree%iwlists, tree%lwlists, tree%iitarget, tree%size, tree%wlists, tree%ntot, &
! 		ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
! 		iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
! 		ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
! 		iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
! 		pot_far_src,  field_far_src,  pot_mpole_src,  field_mpole_src,  pot_near_src,  field_near_src, &
! 		pot_far_targ, field_far_targ, pot_mpole_targ, field_mpole_targ, pot_near_targ, field_near_targ, &
! 		ifnear_only )
    

!     endif ! (solver .eq. SOLVE_DIRECT)

!     ! combine the near and far field for the respective FMM
!     if ( solver .eq. SOLVE_FMM_NEAR .or. solver .eq. SOLVE_FMM_NEAR_FAR) then
!       do i=1, nsources
! 	field_src(1, i) = field_near_src(1,i) + field_far_src(1,i)
! 	field_src(2, i) = field_near_src(2,i) + field_far_src(2,i)
! 	field_src(3, i) = field_near_src(3,i) + field_far_src(3,i)
!       enddo

!       do i=1, ntargets
! 	field_targ(1, i) = field_near_targ(1,i) + field_far_targ(1,i)
! 	field_targ(2, i) = field_near_targ(2,i) + field_far_targ(2,i)
! 	field_targ(3, i) = field_near_targ(3,i) + field_far_targ(3,i)
!       enddo

!     endif

!     ! calculate the accelerations
!     !call calc_accel(accel_src, field_src, charge_src, nsources) 
!     !call calc_accel(accel_targ, field_targ, charge_targ, ntargets)

!     ! store the velocity
!     !f1(1:3, 1:N, 1) 		= y(1:3, 1:N, 2)

!     ! store the acceleartion
!     !f1(1:3, 1:nsources, 2) 	= accel_src(1:3, 1:nsources)
!     !f1(1:3, nsources+1:N, 2) 	= accel_targ(1:3, 1:ntargets)

!     !pdot(1:3, 1:nsources) 	= accel_src(1:3, 1:nsources)
!     !pdot(1:3, nsources+1:N) 	= accel_targ(1:3, 1:ntargets)

!     ! the fancy verlet integrator requires this
!     !qdot = pdot
!     !f1(1:3, 1:N, 2) 	= f1(1:3, 1:N, 1)
    

!     !------------------------------------
!     ! Testing
!     !------------------------------------
!     ! When position is zero. I get NaNs in acceleration.

! !    print *, "pos_src:"
! !    call print_real( y(1:3, 1:nsources, 1), nsources )

! !    print *, "field_src:"
! !    call print_cmplx( field_src, nsources)

! !     print *, "accel_src:"
! !     call print_real( accel_src, nsources )
! !     print *, accel_src(1,1), accel_src(2,1), accel_src(3,1)

!     deallocate( pot_mpole_src )
!     deallocate( pot_mpole_targ )
!     deallocate( field_mpole_src )
!     deallocate( field_mpole_targ )

!     ! the fancy verlet integrator requires this
!     !qdot = pdot

!   end subroutine acceleration


  ! For Testing: my direct solver
  subroutine my_direct(field_src, pot_src, pos_src, nsources, charge_src, field_targ, pot_targ, &
			      pos_targ, ntargets, nthreads)
    
    use omp_lib
    
    implicit none

    ! Function Arguments
    integer(c_int), intent(in), value :: nsources, nthreads, ntargets

    real(pfdp), intent(in) :: pos_src(3, nsources), pos_targ(3, ntargets)

    complex(pfdp), intent(in) ::  charge_src(nsources)

    complex(pfdp), intent(out) :: field_src(3, nsources), pot_src(nsources), field_targ(3, ntargets), &
				   pot_targ(ntargets)

    ! Local variables
    real*8 :: dist, distSq, dx, dy, dz, eMag, ex, ey, ez, qi, qj, pot, K

    integer :: i, j

    !call omp_set_num_threads(nthreads)

    K = 1.0

    !!$omp parallel do &
    !!$omp shared(K, pos_src, charge_src, nsources, pot_src, field_src) &
    !!$omp private(dist, distSq, dx, dy, dz, eMag, ex, ey, ez, qi, qj, pot, i, j) &
    !!$omp& schedule(dynamic)

    do i = 1, nsources
      ! reset the E-field variables
      eMag = 0.0
      ex = 0.0
      ey = 0.0
      ez = 0.0
      pot = 0.0

      qi = Real(charge_src(i))

      do j = 1, nsources
	if (i .ne. j) then
	    ! set the charge for the jth particle
	    qj = Real(charge_src(j))

	    ! calculate the distance
	    dx = pos_src(1,j) - pos_src(1,i)
	    dy = pos_src(2,j) - pos_src(2,i)
	    dz = pos_src(3,j) - pos_src(3,i)
	    distSq = dx*dx + dy*dy + dz*dz
	    dist = sqrt(distSq)

	    ! calculate the x,y, and z E-field
	    eMag = K*qj/distSq
	    ex = ex - eMag*dx/dist	! The minus is because like charges repell NOT attract
	    ey = ey - eMag*dy/dist
	    ez = ez - eMag*dz/dist

	    ! calculate the potential
	    pot = pot + K*qj/dist
	end if
      end do

      ! store the E-field for the ith particle
      field_src(1,i) = ex
      field_src(2,i) = ey
      field_src(3,i) = ez

      ! store the potential for the ith particle
      pot_src(i) = pot
    end do
    !!$omp end parallel do

    !!$omp parallel do &
    !!$omp shared(K, pos_src, pos_targ, charge_src, nsources, ntargets, pot_targ, field_targ) &
    !!$omp private(dist, distSq, dx, dy, dz, eMag, ex, ey, ez, qi, qj, pot, i, j) &
    !!$omp& schedule(dynamic)
    do i = 1, ntargets
      pot  = 0.0
      eMag = 0.0
      ex   = 0.0
      ey   = 0.0
      ez   = 0.0
      qi   = 1.0

      do j = 1, nsources
	! the charge for the jth source
	qj = Real(charge_src(j))

	! calculate the distance
	dx = pos_src(1,j) - pos_targ(1,i)
	dy = pos_src(2,j) - pos_targ(2,i)
	dz = pos_src(3,j) - pos_targ(3,i)
	distSq = dx*dx + dy*dy + dz*dz
	dist = sqrt(distSq)

	! calculate the x,y, and z E-field
	eMag = K*qj/distSq
	ex = ex - eMag*dx/dist	! The minus is because like charges repell NOT attract
	ey = ey - eMag*dy/dist
	ez = ez - eMag*dz/dist

	! calculate the potential
	pot = pot + K*qj/dist
      enddo
      
      ! store the E-field for the ith target
      field_targ(1,i) = ex
      field_targ(2,i) = ey
      field_targ(3,i) = ez

      ! store the potential at the ith targtet
      pot_targ(i) = pot
     
    enddo
    !!$omp end parallel do

  end subroutine my_direct

end module feval
