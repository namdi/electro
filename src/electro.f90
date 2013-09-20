module electro

  use my_type
  use tree

  implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !---------------------------------------------------------
  ! Set solver subroutines
  !---------------------------------------------------------

  !
  ! set solver
  !

  subroutine set_solver(obj)

    type(electro_t), intent(out) :: obj

    integer*4 :: ifsource, iftarget, ifpot, iffld

    ifsource 	= obj%ifsource
    iftarget	= obj%iftarget
    ifpot	= obj%ifpot
    iffld	= obj%iffld

    select case(obj%solver)

    case(SOLVE_DIRECT)
      call set_solver_direct(obj)

    case(SOLVE_FMM)
      call set_solver_direct(obj)

    case(SOLVE_FMM_NEAR_FAR)
      call set_solver_fmm_near_far(obj)

    case(SOLVE_FMM_NEAR)
      call set_solver_fmm_near(obj)

    case default
      stop "ERROR: Unkown solver type."

    end select

  end subroutine set_solver


  !
  ! set flags for the direct solver, which are the same as the full FMM
  !

  subroutine set_solver_direct(obj)

    type(electro_t), intent(out) :: obj

    integer*4 :: ifsource, iftarget, ifpot, iffld

    ifsource 	= obj%ifsource
    iftarget	= obj%iftarget
    ifpot	= obj%ifpot
    iffld	= obj%iffld

    obj%ifpotsrc_far		= 0
    obj%ifpotsrc_near 		= 0

    obj%ifpottarg_far	 	= 0
    obj%ifpottarg_near 		= 0

    obj%iffldsrc_far		= 0
    obj%iffldsrc_near		= 0
      
    obj%iffldtarg_far		= 0
    obj%iffldtarg_near 		= 0

    obj%ifnear_only		= 0

    if (ifpot .eq. 1 .and. ifsource .eq. 1) then
      obj%ifpotsrc = 1
    else
      obj%ifpotsrc = 0
    endif

    if (ifpot .eq. 1 .and. iftarget .eq. 1) then
      obj%ifpottarg = 1
    else
      obj%ifpottarg = 0
    endif

    if (iffld .eq. 1 .and. ifsource .eq. 1) then
      obj%iffldsrc = 1
    else
      obj%iffldsrc = 0
    endif

    if (iffld .eq. 1 .and. iftarget .eq. 1) then
      obj%iffldtarg = 1
    else
      obj%iffldtarg = 0
    endif

  end subroutine set_solver_direct


  !
  ! set flags for the FMM-near-far solver
  !

  subroutine set_solver_fmm_near_far(obj)

    type(electro_t), intent(out) :: obj

    integer*4 :: ifsource, iftarget, ifpot, iffld

    ifsource 	= obj%ifsource
    iftarget	= obj%iftarget
    ifpot	= obj%ifpot
    iffld	= obj%iffld

    obj%ifpotsrc	= 0
    obj%ifpottarg 	= 0

    obj%iffldsrc	= 0
    obj%iffldtarg	= 0
    obj%ifnear_only	= 0

    if (ifpot .eq. 1 .and. ifsource .eq. 1) then
      obj%ifpotsrc_near		= 1
      obj%ifpotsrc_far		= 1
    else
      obj%ifpotsrc_near		= 0
      obj%ifpotsrc_far		= 0
    endif

    if (ifpot .eq. 1 .and. iftarget .eq. 1) then
      obj%ifpottarg_near	= 1
      obj%ifpottarg_far		= 1
    else
      obj%ifpottarg_near	= 0
      obj%ifpottarg_far		= 0
    endif

    if (iffld .eq. 1 .and. ifsource .eq. 1) then
      obj%iffldsrc_near		= 1
      obj%iffldsrc_far		= 1
    else
      obj%iffldsrc_near		= 0
      obj%iffldsrc_far		= 0
    endif

    if (iffld .eq. 1 .and. iftarget .eq. 1) then
      obj%iffldtarg_near 	= 1
      obj%iffldtarg_far 	= 1
    else
      obj%iffldtarg_near 	= 0
      obj%iffldtarg_far 	= 0
    endif

  end subroutine set_solver_fmm_near_far

  !
  ! sets flag for the FMM-near solver
  !

subroutine set_solver_fmm_near(obj)

    type(electro_t), intent(out) :: obj

    integer*4 :: ifsource, iftarget, ifpot, iffld

    ifsource 	= obj%ifsource
    iftarget	= obj%iftarget
    ifpot	= obj%ifpot
    iffld	= obj%iffld

    obj%ifpotsrc	= 0
    obj%ifpotsrc_far	= 0

    obj%ifpottarg 	= 0
    obj%ifpottarg_far	= 0

    obj%iffldsrc	= 0
    obj%iffldsrc_far	= 0

    obj%iffldtarg	= 0
    obj%iffldtarg_far	= 0

    if (ifpot .eq. 1 .and. ifsource .eq. 1) then
      obj%ifpotsrc_near		= 1
    else
      obj%ifpotsrc_near		= 0
    endif

    if (ifpot .eq. 1 .and. iftarget .eq. 1) then
      obj%ifpottarg_near	= 1
    else
      obj%ifpottarg_near	= 0    
    endif

    if (iffld .eq. 1 .and. ifsource .eq. 1) then
      obj%iffldsrc_near		= 1
    else
      obj%iffldsrc_near		= 0     
    endif

    if (iffld .eq. 1 .and. iftarget .eq. 1) then
      obj%iffldtarg_near 	= 1
    else
      obj%iffldtarg_near 	= 0
    endif

  end subroutine set_solver_fmm_near

  !---------------------------------------------------------
  ! Constructor/ Destructor
  !---------------------------------------------------------

  ! 
  ! Construct the electro object
  !
    
  subroutine electro_construct(obj, nsources, ntargets, nthreads, iprec, ifcharge, ifdipole, &
				ifsource, iftarget, solver, ifpot, iffld, lev, charge_src, &
				dipstr_src, dipvec_src, nlevs)
  
    integer, intent(in) :: nsources, ntargets, nthreads, iprec, ifcharge, ifdipole, ifsource, iftarget, & 
			    solver, ifpot, iffld, lev, nlevs
  
    complex*16, target, intent(in) :: charge_src(nsources), dipstr_src(nsources)

    real*8, target, intent(in) :: dipvec_src(3, nsources)

    type(electro_t), intent(out) :: obj
    
    integer :: i
        
    obj%nsources	= nsources
    obj%ntargets	= ntargets
    obj%N		= nsources + ntargets
    obj%nthreads	= nthreads
    obj%shape		= [ 2, obj%N, 3 ]
    obj%size		= size(obj%shape)

    obj%solver		= solver
    obj%ifpot		= ifpot
    obj%iffld		= iffld

    obj%ifcharge	= ifcharge
    obj%ifdipole	= ifdipole

    obj%ifsource	= ifsource
    obj%iftarget	= iftarget

    obj%lev		= lev
    obj%nlevs		= nlevs

    obj%charge_src	=> charge_src
    obj%dipstr_src	=> dipstr_src
    obj%dipvec_src	=> dipvec_src

    ! set the flags for the correct solver
    call set_solver(obj)

    ! tree construction
    if (obj%solver .ne. SOLVE_DIRECT) then
      allocate(obj%tree)
      call tree_construct(obj%tree, nsources, ntargets, iprec)
    endif
    

    allocate( obj%mass_src(nsources) )
    allocate( obj%mass_targ(ntargets) )

    allocate( obj%charge_targ(ntargets) )    
    allocate( obj%pot_src(nsources) )
    allocate( obj%pot_far_src(nsources) )
    allocate( obj%pot_near_src(nsources) )

    allocate( obj%pot_targ(ntargets) )
    allocate( obj%pot_far_targ(ntargets) )
    allocate( obj%pot_near_targ(ntargets) )

    allocate( obj%field_src(3, nsources) )
    allocate( obj%field_far_src(3, nsources) )
    allocate( obj%field_near_src(3, nsources) )

    allocate( obj%field_targ(3, ntargets) )
    allocate( obj%field_far_targ(3, ntargets) )
    allocate( obj%field_near_targ(3, ntargets) )

    allocate( obj%accel_src(3, nsources) )
    allocate( obj%accel_targ(3, ntargets) )

    ! set the masses and target charges
    do i=1,nsources
      obj%mass_src(i) = 1.0
    enddo

    do i=1,ntargets
      obj%mass_targ(i) = 1.0
      obj%charge_targ(i) = cmplx(1.0)
    enddo

  end subroutine electro_construct

  !
  ! Destruct electro object
  !

  subroutine electro_destruct(obj)

    type(electro_t), intent(out) :: obj

    ! deallocate tree
    if (obj%solver .ne. SOLVE_DIRECT) then
      call tree_destruct(obj%tree)
    endif

    deallocate( obj%mass_src )
    deallocate( obj%mass_targ)
    deallocate( obj%charge_targ )
    
    
    deallocate( obj%pot_src )
    deallocate( obj%pot_far_src )
    deallocate( obj%pot_near_src )

    deallocate( obj%pot_targ )
    deallocate( obj%pot_far_targ )
    deallocate( obj%pot_near_targ )

    deallocate( obj%field_src )
    deallocate( obj%field_far_src )
    deallocate( obj%field_near_src )

    deallocate( obj%field_targ )
    deallocate( obj%field_far_targ )
    deallocate( obj%field_near_targ)

    deallocate( obj%accel_src )
    deallocate( obj%accel_targ )

    nullify( obj%mass_src )
    nullify( obj%mass_targ )
    
    nullify( obj%charge_src )
    nullify( obj%charge_targ )
    nullify( obj%dipvec_src )
    nullify( obj%dipstr_src )

    nullify( obj%pot_src )
    nullify( obj%pot_far_src )
    nullify( obj%pot_near_src )

    nullify( obj%pot_targ )
    nullify( obj%pot_far_targ )
    nullify( obj%pot_near_targ )

    nullify( obj%field_src )
    nullify( obj%field_far_src )
    nullify( obj%field_near_src )
    
    nullify( obj%field_targ )
    nullify( obj%field_far_targ )
    nullify( obj%field_near_targ )

    nullify( obj%accel_src )
    nullify( obj%accel_targ )

  end subroutine electro_destruct

end module electro