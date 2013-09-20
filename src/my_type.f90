module my_type
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none

  integer, parameter :: LADDER_MAX = 200

  integer, parameter :: SOLVE_FMM_NEAR 		= 1
  integer, parameter :: SOLVE_FMM_NEAR_FAR	= 2
  integer, parameter :: SOLVE_FMM   		= 3
  integer, parameter :: SOLVE_DIRECT 		= 4
  !
  ! Tree object
  !	size: size of the box on level 0

  type :: tree_t

    integer*4 :: ladder(2, LADDER_MAX), box(20)

    real*8 :: center_temp(3), corners_temp(3,8), epsfmm

    integer*4 :: nsources, ntargets, N, iisource, iitarget, iwlists, nlev, &
	      nboxes, nbox, ntot, iprec, lwlists, size, lused7
  
    real*8, pointer :: wlists(:)

    integer*4, pointer :: isource_tree(:), itarget_tree(:), boxes(:,:), &
			  centers(:,:), corners(:,:,:)

  end type tree_t

  !
  ! Electro object
  !

  type :: electro_t

    integer*4		:: nsources, ntargets, N, shape(3), size, nthreads
  
    ! iffld whether or not any of the fields will be calculated
    integer*4		:: solver, ifpot, iffld, ifsource, iftarget, ifcharge, ifdipole, ifnear_only, lev, nlevs

    integer*4		:: ifpotsrc, ifpotsrc_far, ifpotsrc_near, ifpottarg, ifpottarg_far, ifpottarg_near

    integer*4		:: iffldsrc, iffldsrc_far, iffldsrc_near, iffldtarg, iffldtarg_far, iffldtarg_near

    integer*4		:: dofmm, dofmm_multi, dodirect, dofmm_near

    type(tree_t), pointer	:: tree
    
    real(pfdp), pointer		:: mass_src(:), mass_targ(:), dipvec_src(:,:), accel_src(:,:), accel_targ(:,:)

    complex(pfdp), pointer	:: pot_src(:), pot_far_src(:), pot_near_src(:), dipstr_src(:), charge_src(:), &
				    pot_targ(:), pot_far_targ(:), pot_near_targ(:), charge_targ(:)

    complex(pfdp), pointer 	:: field_src(:,:), field_near_src(:,:), field_far_src(:,:), field_targ(:,:), &
				    field_far_targ(:,:), field_near_targ(:,:)

  end type electro_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module my_type
