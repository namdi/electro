! This file contains some fortran subroutines to be used in python for the fmm
! these functions are from lfmm3dpart.f

module fortfmm
  use iso_c_binding 	! makes the fortran code look like c
contains
  
  subroutine create_tree_wrap(iprec, source, nsource, target, ntarget, & 
		center, laddr, idata, idata_len, &
		ddata, ddata_len, wlists, ntot) bind(c, name = "create_tree_wrap")
    
    implicit none
    ! idata is an array that contains integer data so that I can return it to the python functions.
    ! it is a hack to get around the fact that python does not allow scalars to be returned
    ! by reference.
    ! idata(1) = ier	idata(2) = nlev 	idata(3) = nboxes	idata(4) = nbox		
    ! idata(5) = lused7	idata(6) = iisource	idata(7) = iwlists	idata(8) = lwlists
    ! idata(9) = iitarget 	idata(10) = not

    ! ddata(1) = epsfmm		ddata(2) = size

    integer(c_int), intent(in), value :: iprec, nsource, ntarget, idata_len, ddata_len, ntot

    real(c_double), intent(in) :: source(3, nsource), target(3, ntarget), center(3)

    integer(c_int), intent(out) :: laddr(2, 200), idata(idata_len), wlists(ntot)

    !integer(c_int) :: ier, nlev, nboxes, nbox, lused7
  
    real(c_double), intent(out) :: ddata(ddata_len)

    real(c_double) :: epsfmm, size

    integer(c_int) :: ier, nlev, nboxes, nbox, lused7, iisource, iwlists, lwlists, iitarget

    !real(c_double), allocatable :: wlists(:)

    !write(*,*) "In create_tree_wrap() fortfmm.f90"
    call create_tree(ier, iprec, source, nsource, target, ntarget, center, laddr, nlev, nboxes, &
		      nbox, epsfmm, lused7, iisource, iwlists, lwlists, iitarget, size, wlists, ntot)

    ! store integer data for output
    idata(1) = ier
    idata(2) = nlev
    idata(3) = nboxes
    idata(4) = nbox
    idata(5) = lused7
    idata(6) = iisource
    idata(7) = iwlists
    idata(8) = lwlists
    idata(9) = iitarget

    ddata(1) = epsfmm
    ddata(2) = size

    !write(*,*) "in fortfmm.f90 create_tree_wrap() epsfmm ", epsfmm
    
  end subroutine create_tree_wrap

  subroutine pre_create_tree_wrap(iprec, source, nsource, target, ntarget, & 
		idata, idata_len) bind(c, name = "pre_create_tree_wrap")
    
    implicit none
    ! idata is an array that contains integer data so that I can return it to the python functions.
    ! it is a hack to get around the fact that python does not allow scalars to be returned
    ! by reference.
    ! idata(1) = ier	idata(2) = tot

    integer(c_int), intent(in), value :: iprec, nsource, ntarget, idata_len

    real(c_double), intent(in) :: source(3, nsource), target(3, ntarget)

    integer(c_int), intent(out) :: idata(idata_len)

    integer(c_int) :: ier, nlev, nboxes, nbox, lused7, iisource, iwlists, lwlists, iitarget, ntot

    !write(*,*) "In pre_create_tree_wrap() fortfmm.f90"

    ! the sole purpose of this function is to create a tree-structure so that I can find the 
    ! length of the tree-structure which is "ntot"
    call pre_create_tree(ier, iprec, source, nsource, target, ntarget, nlev, nbox, ntot)

    ! store integer data for output
    idata(1) = ier
    idata(2) = ntot

  end subroutine pre_create_tree_wrap


  subroutine get_box_info_wrap(ibox, box, center, corners, w, idata, idata_len, iwlists) &
	      bind(c, name = "get_box_info_wrap")

    implicit none

    integer(c_int), intent(in), value :: ibox, idata_len, iwlists

    real(c_double), intent(in) :: w(*)

    integer(c_int), intent(out) :: box(20), idata(idata_len)

    real(c_double), intent(out) :: center(3), corners(3,8)

    integer(c_int) :: ier

    !idata(1) = ier

    ier = idata(1)

    call get_box_info(ier, ibox, box, center, corners, w(iwlists))

    idata(1) = ier

  end subroutine get_box_info_wrap

  subroutine get_winfo_wrap(my_wlists, wlists, idx, n) bind(c, name="get_winfo_wrap")

    implicit none

    integer(c_int), intent(in), value :: idx, n

    real(c_double), intent(in) :: wlists(*)

    integer(c_int), intent(out) :: my_wlists(*)

    call get_winfo(my_wlists, wlists(idx), n)
  
  end subroutine get_winfo_wrap

  subroutine fmm_wrap(iprec, nsource, source, ifcharge, charge, ifdipole, &
		      dipstr, dipvec, ifpot, pot, iffld, fld, &
		      ntarget, target, ifpottarg, pottarg, iffldtarg, &
		      fldtarg, laddr, idata, idata_len, ddata, ddata_len, &
		      wlists, ntot, &
		      ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
		      iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
		      ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
		      iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
		      potsrc_far,  fldsrc_far, potsrc_mpole, fldsrc_mpole, &
		      potsrc_near, fldsrc_near, &
		      pottarg_far, fldtarg_far, &
		      pottarg_mpole, fldtarg_mpole, pottarg_near, &
		      fldtarg_near, ifnear_only) bind(c, name = "fmm_wrap")

    ! we should copy the data from the real and imaginary parts in parallel
    ! complex(c_double_complex)

    use omp_lib !use this for openmp functionality

    implicit none

    ! Inputs into fmm_wrap
    integer(c_int), intent(in), value :: iprec, ifcharge, ifdipole, ifpot, iffld, idata_len, &
					  nsource, ntarget, ifpottarg, iffldtarg, ddata_len, &
					  ntot, &
					  ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
					  iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
					  ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
					  iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
					  ifnear_only
					  
    integer(c_int), intent(in) :: laddr(2, 200)

    real(c_double), intent(in) :: source(3, nsource), dipvec(3, nsource), target(3, nsource)

    real(c_double), intent(in) :: ddata(ddata_len)
    
    real(c_double_complex), intent(in) :: charge(nsource), dipstr(nsource) 

    ! ouputs

    real(c_double_complex), intent(out) :: pot(nsource), fld(3, nsource), pottarg(ntarget), &
					   fldtarg(3, ntarget), & 
					   potsrc_far(nsource),   fldsrc_far(3, nsource), &
					   potsrc_mpole(nsource), fldsrc_mpole(3, nsource), &
					   potsrc_near(nsource),  fldsrc_near(3, nsource), &
					   pottarg_far(ntarget), &
					   fldtarg_far(3, ntarget), pottarg_mpole(ntarget), &
					   fldtarg_mpole(3, ntarget), &
					   pottarg_near(ntarget), &
					   fldtarg_near(3, ntarget)
    
    integer(c_int), intent(out) ::  idata(idata_len)
    
    real(c_double), intent(out) :: wlists(ntot)
    !integer(c_int), intent(out) :: wlists(ntot)

    ! arguments
	
    integer(c_int) :: ier, nlev, nboxes, nbox, lused7, nthread, i, j, iisource, iwlists, lwlists, iitarget

    real(c_double) :: epsfmm, size

    
    ! idata(1) = ier	idata(2) = nlev	idata(3) = nboxes	idata(4) = nbox
    ! idata(5) = lused7	idata(6) = nthread	idata(7) = iisource	idata(8) = iwlists
    ! idata(9) = lwlists	idata(10) = iitarget

    ! ddata(1) = epsfmm		ddata(2) = size

    ier	 	= idata(1)
    nlev 	= idata(2)
    nboxes 	= idata(3)
    nbox 	= idata(4)
    lused7 	= idata(5)
    nthread	= idata(6)
    iisource	= idata(7)
    iwlists	= idata(8)
    lwlists	= idata(9)
    iitarget	= idata(10)

    epsfmm	= ddata(1)
    size	= ddata(2)
    
    !write(*,*) "In fortfmm.f90"
    !write(*,*) "nthread: " , nthread, "iprec: ", iprec
    !write(*,*) "ifcharge: ", ifcharge, "ifdipole: ", ifdipole
    !write(*,*) "ifpot: " , ifpot, "iffld: ", iffld
    !write(*,*) "idata_len: " , idata_len, "nsource: " , nsource
    !write(*,*) "ntarget: " , ntarget, "epsfmm: ", epsfmm
    !write(*,*) "iisource: " , iisource, "iitarget: ", iitarget
    !write(*,*) "iwlists: " , iwlists, "lwlists: ", lwlists
    !write(*,*) "ifnear_only(fortfmm.f90): ", ifnear_only
    !-------------------------------------------------------
    ! run the fmm calculation
    !-------------------------------------------------------
    !write(*,*) "About to enter fmm() "

    call fmm(ier, iprec, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, &
		ifpot, pot, iffld, fld, ntarget, target, ifpottarg, pottarg, iffldtarg, &
		fldtarg, laddr, nlev, nboxes, nbox, epsfmm, lused7, nthread, iisource, &
		iwlists, lwlists, iitarget, size, wlists, ntot, &
		ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
		iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
		ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
		iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
		potsrc_far, fldsrc_far, potsrc_mpole, fldsrc_mpole, potsrc_near, fldsrc_near, &
		pottarg_far, fldtarg_far, &
		pottarg_mpole, fldtarg_mpole, pottarg_near, fldtarg_near, &
		ifnear_only)

    !write(*,*) "Just finished fmm() in fortfmm.f90"
    
    ! For integer output
    idata(1)	= ier
    idata(2)	= nlev
    idata(3) 	= nboxes
    idata(4)	= nbox
    idata(5)	= lused7
    idata(6)	= nthread
    
  !write(*,*) "leaving fmm_wrap in fortfmm.f90"
  end subroutine fmm_wrap


  !subroutine direct_wrap(nsources, source, ifcharge, charge, ifdipole, dipstr, dipvec, & 
			  !ifpot, pot, iffld, fld, ntargets, target, ifpottarg, pottarg, &
			  !iffldtarg, fldtarg, nthreads) bind(c, name="direct_wrap")

  subroutine direct_wrap( idata, idata_len, source, charge, dipstr, dipvec, & 
			  pot, fld, target, pottarg, fldtarg ) &
			  bind(c, name="direct_wrap")

  use omp_lib

  implicit none

  !idata(1) = ifcharge	idata(2) = ifdipole	idata(3) = ifpot	idata(4) = iffld
  !idata(5) = ifpottarg	idata(6) = iffldtarg	idata(7) = nthreads	idata(8) = nsources
  !idata(9) = ntargets

  !integer(c_int), intent(in), value :: nsources, ifcharge, ifdipole, ifpot, iffld, ntargets, &
					!ifpottarg, iffldtarg, nthreads

  integer(c_int), intent(in), value :: idata_len

  real(c_double_complex), intent(in) :: charge(*), dipstr(*)
  
  integer(c_int), intent(in) :: idata(*)

  real(c_double), intent(in) :: source(3, *), target(3, *), dipvec(3, *)
  
  real(c_double_complex), intent(out) :: pot(*), fld(3, *), pottarg(*), &
					  fldtarg(3, *)

  integer(c_int) :: nsources, ifcharge, ifdipole, ifpot, iffld, ntargets, ifpottarg, iffldtarg, nthreads

  ifcharge 	= idata(1)
  ifdipole 	= idata(2)
  ifpot		= idata(3)
  iffld		= idata(4)
  ifpottarg	= idata(5)
  iffldtarg	= idata(6)
  nthreads	= idata(7)
  nsources	= idata(8)
  ntargets	= idata(9)

  !write(*,*) "In direct_wrap() of fortfmm.f90"
  !write(*,*) "ifcharge, ", ifcharge, "ifdipole, ", ifdipole
  !write(*,*) "ifpot, ", ifpot, "iffld, ", iffld
  !write(*,*) "ifpottarg, ", ifpottarg, "iffldtarg, ", iffldtarg
  !write(*,*) "nthreads, ", nthreads, "nsources, ", nsources
  !write(*,*) "ntargets, ", ntargets

  call direct( nsources, source, ifcharge,charge, ifdipole, dipstr, dipvec, ifpot, pot, & 
	      iffld, fld, ntargets, target, ifpottarg, pottarg, iffldtarg, fldtarg, nthreads )

  end subroutine direct_wrap

end module fortfmm

