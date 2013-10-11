  !
  ! Create the tree
  !

  subroutine create_tree(ier, iprec, source, nsource, target, ntarget, center, laddr, nlev, nboxes, &
			nbox, epsfmm, lused7, iisource, iwlists, lwlists, iitarget, size, wlists, ntot)
    ! This method is a snipet of lfmm3dparttarg() in lfmm3dpart.f
    implicit real *8 (a-h,o-z)

    ! function arguments
    integer, intent(in) :: iprec, nsource, ntarget, ntot

    real *8, intent(in) :: source(3, nsource), target(3, ntarget)

    integer, intent(out) :: laddr(2, 200), center(3) !, wlists(ntot)

    integer, intent(out) :: ier, nlev, lused7, nboxes, nbox, iisource, iwlists, lwlists, iitarget

    real *8, intent(out) :: epsfmm, size, wlists(ntot)

    ! method variables
    !integer :: ntot
    !real *8, allocatable :: wlists(:)
	
    ! nbox - the maximum number of points in a box on the finest level     

	ier=0 
	lused7 = 0

	iisource = 0
	iitarget = 0
	iwlists  = 0
	lwlists  = 0
	nboxes 	 = 0
	center(1) = 0.0
	center(2) = 0.0
	center(3) = 0.0

	! Namdi: created the following subroutines
	call set_fmmtolerance(iprec, epsfmm)
	
	call set_sources_per_box(iprec, nsource, ntarget, nbox)

        !ntot = 100*(nsource+ntarget)+10000
	
	
	! Namdi: IMPORTANT:
	! since I use pre_create_tree() and already know the length (ntot),
	! THIS SHOULD NOT BE A LOOP

	! try to allocate wlists at most 10 times
        !do ii = 1,10
           !allocate (wlists(ntot))
	   !write(*,*) "ii: ", ii
            call lfmm3dparttree(ier,iprec, &
             nsource,source,ntarget,target, &
             nbox,epsfmm,iisource,iitarget,iwlists,lwlists, &
             nboxes,laddr,nlev,center,size, &
             wlists,ntot,lused7)
           !if (ier.ne.0) then
              !deallocate(wlists)
              !ntot = ntot*1.5
              !call prinf(' increasing allocation, ntot is *',ntot,1)
           !else
             !goto 1200
           !endif
        !enddo
!1200    continue
	  


        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return
        endif

    !write(*,*) "iwlists: ", iwlists, "lwlists: ", lwlists
    !write(*,*) "size: " , size

  end subroutine create_tree


  !    
  ! this function finds ntot, the length of wlists when we actually create the tree
  !
  subroutine pre_create_tree(ier, iprec, source, nsource, target, ntarget, nlev, &
			nbox, ntot)
    ! This method is a snipet of lfmm3dparttarg() in lfmm3dpart.f
    implicit real *8 (a-h,o-z)

    ! function arguments
    integer, intent(in) :: iprec, nsource, ntarget

    real *8, intent(in) :: source(3, nsource), target(3, ntarget)

    integer, intent(out) :: ier, nlev, nbox, ntot

    ! method variables
    real *8, allocatable :: wlists(:)

    real *8 :: size, epsfmm, center(3)

    integer :: lwlists, iwlists, iisource, iitarget, lused7, nboxes

    integer :: laddr(2, 200)
	
    ! nbox - the maximum number of points in a box on the finest level 
 
	ier=0 
	lused7 = 0

	iisource = 0
	iitarget = 0
	iwlists  = 0
	lwlists  = 0
	nboxes 	 = 0
	center(1) = 0.0
	center(2) = 0.0
	center(3) = 0.0

	! Namdi: created the following subroutines
	call set_fmmtolerance(iprec, epsfmm)
	call set_sources_per_box(iprec, nsource, ntarget, nbox)

        ntot = 100*(nsource+ntarget)+10000

	! try to allocate wlists at most 10 times
        do ii = 1,10
           allocate (wlists(ntot))
	   !write(*,*) "ii: ", ii
            call lfmm3dparttree(ier,iprec, &
             nsource,source,ntarget,target, &
             nbox,epsfmm,iisource,iitarget,iwlists,lwlists, &
             nboxes,laddr,nlev,center,size, &
             wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
             goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return
        else
	  deallocate(wlists)	! Namdi: added the else clause
	endif

  end subroutine pre_create_tree


  !
  ! this function sets the fmm tolerance
  !
  subroutine set_fmmtolerance(iprec, epsfmm) 

    implicit none

    integer, intent(in) :: iprec

    real*8, intent(out) :: epsfmm
  
    !set fmm tolerance based on iprec flag.
    if( iprec .eq. -2 ) epsfmm=.5d-0 
    if( iprec .eq. -1 ) epsfmm=.5d-1
    if( iprec .eq. 0 ) epsfmm=.5d-2
    if( iprec .eq. 1 ) epsfmm=.5d-3
    if( iprec .eq. 2 ) epsfmm=.5d-6
    if( iprec .eq. 3 ) epsfmm=.5d-9
    if( iprec .eq. 4 ) epsfmm=.5d-12
    if( iprec .eq. 5 ) epsfmm=.5d-15
    if( iprec .eq. 6 ) epsfmm=0
      
    !if (ifprint .ge. 1) call prin2('epsfmm=*',epsfmm,1)

  end subroutine set_fmmtolerance

  !
  ! this sets the maximum number of sources per box 
  !

  subroutine set_sources_per_box(iprec, nsource, ntarget, nbox)
  !set criterion for box subdivision (number of sources per box)

    implicit none

    integer, intent(in) :: iprec, nsource, ntarget
  
    integer, intent(out) :: nbox

    if( iprec .eq. -2 ) nbox=40*1.0
    if( iprec .eq. -1 ) nbox=50*1.0
    if( iprec .eq. 0 ) nbox=80*1.0
    if( iprec .eq. 1 ) nbox=160*1.0
    if( iprec .eq. 2 ) nbox=400*1.0
    if( iprec .eq. 3 ) nbox=800*1.0 ! should be 800
    if( iprec .eq. 4 ) nbox=1200*1.0
    if( iprec .eq. 5 ) nbox=1400*1.0
    if( iprec .eq. 6 ) nbox=nsource+ntarget

    !if (ifprint .ge. 1) call prinf('nbox=*',nbox,1)

  end subroutine set_sources_per_box

  !
  ! This actually runs the FMM
  !

  subroutine fmm(ier, iprec, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, &
		ifpot, pot, iffld, fld, ntarget, target, ifpottarg, pottarg, iffldtarg, &
		fldtarg, laddr, nlev, nboxes, nbox, epsfmm, lused7, nthread, iisource, &
		iwlists, lwlists, iitarget, size, wlists, ntot, &
		ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
		iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
		ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
		iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
		potsrc_far,  fldsrc_far,  potsrc_mpole,  fldsrc_mpole,  potsrc_near,  fldsrc_near, &
		pottarg_far, fldtarg_far, pottarg_mpole, fldtarg_mpole, pottarg_near, fldtarg_near, &
		ifnear_only )

  ! This algorithm runs the fmm method
  ! This program is based on lfmm3dparttarg() of lfmm3dpart.f

!C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       
!       Laplace FMM in R^3: evaluate all pairwise particle
!       interactions (ignoring self-interaction) 
!       and interactions with targets.
!
!       We use (1/r) for the Green's function,
!       without the (1/4 pi) scaling.  Self-interactions are not included.
!   
!       This is primarily a memory management code. 
!       The actual work is carried out in subroutine lfmm3dparttargmain.
!
!       INPUT PARAMETERS:
!
!       iprec:     FMM precision flag
!
!                 -2 => tolerance =.5d0
!                 -1 => tolerance =.5d-1
!                  0 => tolerance =.5d-2
!                  1 => tolerance =.5d-3
!                  2 => tolerance =.5d-6
!                  3 => tolerance =.5d-9
!                  4 => tolerance =.5d-12
!                  5 => tolerance =.5d-15
!
!       nsource          : number of sources                           (integer)
!       source(3,nsource): source locations 			       (real *8)  
!       ifcharge         : charge computation flag                     (integer)
!                          ifcharge = 1 => include charge contribution
!                                     otherwise do not
!       charge(nsource)  : charge strengths                        (complex *16)
!       ifdipole         : dipole computation flag                     (integer)
!                          ifdipole = 1 =>  include dipole contribution
!                                     otherwise do not
!       dipstr(nsource)  : dipole strengths                        (complex *16) 
!       dipvec(3,nsource): dipole orientation vectors                  (real *8) 
!
!       ifpot            : potential flag                              (integer)
!                          (1=compute potential, otherwise no)
!       iffld            : field flag                                  (integer) 
!                          (1=compute field, otherwise no)
!       ntarget          : number of targets                           (integer)  
!       target(3,ntarget): target locations                            (real *8) 
!       ifpottarg        : target potential flag                       (integer)
!                          (1=compute potential, otherwise no)
!       iffldtarg        : target field flag                           (integer)
!                          (1=compute field, otherwise no)
!
!       OUTPUT PARAMETERS:
!
!       ier   =  error return code
!                ier=0     =>  normal execution
!                ier=4     =>  cannot allocate tree workspace
!                ier=8     =>  cannot allocate bulk FMM  workspace
!                ier=16    =>  cannot allocate mpole expansion
!                              workspace in FMM
!
!       pot(nsource)       : potential at source locations         (complex *16)
!       fld(3,nsource)     : field (-gradient) at source locations (complex *16)
!       pottarg(ntarget)   : potential at target locations         (complex *16)
!       fldtarg(3,ntarget) : field (-gradient) at target locations (complex *16) 
!-----------------------------------------------------------------------
!
	! commenting out the cf2py codes
!cf2py   intent(out) ier
!cf2py   intent(in) iprec
!cf2py   intent(in) nsource, source
!cf2py   intent(in) ifcharge,charge
!cf2py   check(!ifcharge || (shape(charge,0) == nsource))  charge
!cf2py   depend(nsource)  charge
!cf2py   intent(in) ifdipole,dipvec,dipstr
!cf2py   check(!ifdipole || (shape(dipstr,0) == nsource))  dipstr
!cf2py   depend(nsource)  dipstr
!cf2py   intent(in) ifpot,iffld
!cf2py   intent(out) pot,fld
!cf2py   intent(in) ifpottarg, iffldtarg
!cf2py   intent(in) target
!cf2py   intent(in) ntarget
!cf2py   check((!ifpottarg && !iffldtarg) || (shape(target,0)==3 && shape(target,1) == ntarget))  target
!cf2py   check((!ifpottarg) || (shape(pottarg,0)==ntarget))  pottarg
!cf2py   check((!iffldtarg) || (shape(fldtarg,0)==3 && shape(fldtarg,1) == ntarget))  fldtarg
!
!       (F2PY workaround: pottarg, fldtarg must be input because f2py
!       refuses to allocate zero-size output arrays.)
!
!cf2py   intent(in,out) pottarg,fldtarg
!
	! Namdi: added the following
	use omp_lib

        implicit real *8 (a-h,o-z)	
		      
	!-------------------------------------
	! Namdi: The follwing is code that I've written
	!------------------------------------
	integer, intent(in) :: iprec, ifcharge, ifdipole, ifpot, iffld, ntarget, ntot, &
				ifpottarg, iffldtarg, nlev, nboxes, nbox, nthread, &
				ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
				iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
				ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
				iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
			        ifnear_only

	integer, intent(in) :: laddr(2, 200)

	real *8, intent(in) :: source(3, nsource), dipvec(3, nsource), target(3, nsource), epsfmm, size
    
	complex *16, intent(in) :: charge(nsource), dipstr(nsource)
	
	integer, intent(out) :: ier, lused7, iisource, iwlists, lwlists, iitarget

	complex *16, intent(out) :: pot(nsource), fld(3, nsource), pottarg(ntarget), fldtarg(3, ntarget), &
				    potsrc_far(nsource), fldsrc_far(3, nsource), &
				    potsrc_mpole(nsource), fldsrc_mpole(3, nsource), &
				    potsrc_near(nsource), fldsrc_near(3, nsource), &

				    pottarg_far(ntarget), fldtarg_far(3, ntarget), &
				    pottarg_mpole(ntarget), fldtarg_mpole(3, ntarget), &
				    pottarg_near(3, ntarget), fldtarg_near(3, ntarget)

	real *8, intent(out) :: wlists(ntot)
	!integer, intent(out) :: wlists(ntot)

        !dimension source(3,nsource)
        !complex *16 charge(nsource)
        !complex *16 dipstr(nsource)
        !dimension dipvec(3,nsource)
        complex *16 :: ima
        !complex *16 pot(nsource)
        !complex *16 fld(3,nsource)
        !dimension target(3,nsource)
        !complex *16 pottarg(ntarget)
        !complex *16 fldtarg(3,ntarget)
!
        dimension timeinfo(10)
!
!     Note: various arrays dimensioned here to 200.
!     That allows for 200 evels of refinment, which is 
!     more than enough for any non-pathological case.
!
 
        !dimension laddr(2,200)		passed into this function
        dimension bsize(0:200)
        dimension nterms(0:200)
        integer box(20)
        integer box1(20)
        dimension scale(0:200)
        dimension center(3)
        dimension center0(3),corners0(3,8)
        dimension center1(3),corners1(3,8)
        real *8, allocatable :: w(:)
        !real *8, allocatable :: wlists(:) Already declared up top
        real *8, allocatable :: wrmlexp(:)
        complex *16 :: ptemp,ftemp(3)
       
        data ima/(0.0d0,1.0d0)/
       
	!-----------------------------------
	! Namdi
	!write(*,*) "In fmm() in fortfmm_utils.f90"


        ier=0
        lused7 = 0
       
        done=1
        pi=4*atan(done)

!     ifprint is an internal information printing flag. 
!     Suppressed if ifprint=0.
!     Prints timing breakdown and other things if ifprint=1.
       
!        ifprint=1
	ifprint =0
      
!     set fmm tolerance based on iprec flag.

       !if( iprec .eq. -2 ) epsfmm=.5d-0 
        !if( iprec .eq. -1 ) epsfmm=.5d-1
        !if( iprec .eq. 0 ) epsfmm=.5d-2
        !if( iprec .eq. 1 ) epsfmm=.5d-3
        !if( iprec .eq. 2 ) epsfmm=.5d-6
        !if( iprec .eq. 3 ) epsfmm=.5d-9
        !if( iprec .eq. 4 ) epsfmm=.5d-12
        !if( iprec .eq. 5 ) epsfmm=.5d-15
        !if( iprec .eq. 6 ) epsfmm=0
       
        if (ifprint .ge. 1) call prin2('epsfmm=*',epsfmm,1)
	

!     set criterion for box subdivision (number of sources per box)

        !if( iprec .eq. -2 ) nbox=40*1.0
        !if( iprec .eq. -1 ) nbox=50*1.0
        !if( iprec .eq. 0 ) nbox=80*1.0
        !if( iprec .eq. 1 ) nbox=160*1.0
        !if( iprec .eq. 2 ) nbox=400*1.0
        !if( iprec .eq. 3 ) nbox=800*1.0
        !if( iprec .eq. 4 ) nbox=1200*1.0
        !if( iprec .eq. 5 ) nbox=1400*1.0
        !if( iprec .eq. 6 ) nbox=nsource+ntarget

        if (ifprint .ge. 1) call prinf('nbox=*',nbox,1)


!     create oct-tree data structure

        !t1=second()
!$        t1=omp_get_wtime()
        !ntot = 100*(nsource+ntarget)+10000
        !do ii = 1,10
           !allocate (wlists(ntot))
           !call lfmm3dparttree(ier,iprec,
     !$        !nsource,source,ntarget,target,
     !$        !nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     !$        !nboxes,laddr,nlev,center,size,
     !$        !wlists,ntot,lused7)
           !if (ier.ne.0) then
              !deallocate(wlists)
              !ntot = ntot*1.5
              !call prinf(' increasing allocation, ntot is *',ntot,1)
           !else
             !goto 1200
           !endif
        !enddo
!1200    !continue
        !if (ier.ne.0) then
           !call prinf(' exceeded max allocation, ntot is *',ntot,1)
           !ier = 4
           !return
        !endif
        !t2=second()
!$        t2=omp_get_wtime()
        if( ifprint .eq. 1 ) call prin2('time in d3tstrcr=*',t2-t1,1)

!     lused7 is counter that steps through workspace,
!     keeping track of total memory used.
        lused7=1
        do i = 0,nlev
        scale(i) = 1.0d0
        enddo
       
        if (ifprint .ge. 1) call prin2('scale=*',scale,nlev+1)
       

!       carve up workspace further

!     isourcesort is pointer for sorted source coordinates
!     itargetsort is pointer for sorted target locations
!     ichargesort is pointer for sorted charge densities
!     idipvecsort is pointer for sorted dipole orientation vectors
!     idipstrsort is pointer for sorted dipole densities

        isourcesort = lused7 + 5
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 3*nsource
          ldipstr = 2*nsource
        else
          ldipvec = 3
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr

!       ... allocate the potential and field arrays

        ipot = lused7
        lpot = 2*nsource
        lused7=lused7+lpot
       
        ifld = lused7
        if( iffld .eq. 1) then
        lfld = 2*(3*nsource)
        else
        lfld=6
        endif
        lused7=lused7+lfld
      
        ipottarg = lused7
        lpottarg = 2*ntarget
        lused7=lused7+lpottarg
       
        ifldtarg = lused7
        if( iffldtarg .eq. 1) then
        lfldtarg = 2*(3*ntarget)
        else
        lfldtarg=6
        endif
        lused7=lused7+lfldtarg
      
        if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)

!       based on FMM tolerance, compute expansion lengths nterms(i)
!	------------------------------------------------------------------
!	Namdi: Control the number of terms in the expansions
!	--------------------------------------------------------------------      
        nmax = 0
        call l3dterms(epsfmm, nterms_lap, ier)
        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           nterms(i)=nterms_lap
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo

        if (ifprint .ge. 1) call prinf('nterms=*',nterms,nlev+1)
        if (ifprint .ge. 1) call prinf('nmax=*',nmax,1)

!     Multipole and local expansions will be held in workspace
!     in locations pointed to by array iaddr(2,nboxes).

!     iiaddr is pointer to iaddr array, itself contained in workspace.
!     imptemp is pointer for single expansion (dimensioned by nmax)
   
!       ... allocate iaddr and temporary arrays

        iiaddr = lused7 
        imptemp = iiaddr + 2*nboxes
        lmptemp = (nmax+1)*(2*nmax+1)*2
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace  lused7 is *', &
		      lused7,1)
           ier = 8
           return
        endif

!     reorder sources, targets so that each box holds
!     contiguous list of source/target numbers.

        call l3dreorder(nsource,source,ifcharge,charge,wlists(iisource), &
			ifdipole,dipstr,dipvec, &
			w(isourcesort),w(ichargesort),w(idipvecsort),w(idipstrsort)) 

        call l3dreordertarg(ntarget,target,wlists(iitarget), &
			    w(itargetsort))

        if (ifprint .ge. 1) call prinf('finished reordering=*',ier,1)
        if (ifprint .ge. 1) call prinf('ier=*',ier,1)
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('nlev=*',nlev,1)
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('lused7=*',lused7,1)

!     allocate memory need by multipole, local expansions at all
!     levels
!     irmlexp is pointer for workspace need by various fmm routines,

        call l3dmpalloc(wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)

        if (ifprint .ge. 1) call prinf(' lmptot is *',lmptot,1)
       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace lused7 is *', &
		      lused7,1)
           ier = 16
           return
        endif

	 ! The following is commented out by the original source code
!        do i=lused7+1,lused7+1+100
!        w(i)=777
!        enddo

!     Memory allocation is complete. 

        t1=second()
!$        t1=omp_get_wtime()

	! Namdi: note that I am calling my_lfmm3dparttargmain.
	! my_lfmm3dparttargmain is declared in the source folder to avoid
	! .f to .f90 conversion problems
	
	!write(*,*) "ifnear_only (fortfmm_utils) ", ifnear_only

!--------------------------------------------------------------------------
!     Call main fmm routine. There are, unfortunately, a lot
!     of parameters here. ifevalfar and ifevalloc determine
!     whether far field and local fields (respectively) are to 
!     be evaluated. Setting both to 1 means that both will be
!     computed (which is the normal scenario).
!--------------------------------------------------------------------------
        ifevalfar=1
        ifevalloc=1
	
	if (ifnear_only .eq. 1) then
	  ifevalfar= 0
	  ifevalloc= 1
	endif
!--------------------------------------------------------------------------

	!write(*,*)
	!write(*,*) "Before the calc, in fortfmm_utils.f90 fld(1,1): ", fld(3,1)

	call my_lfmm3dparttargmain(ier,iprec, &
          ifevalfar,ifevalloc, &
          nsource,w(isourcesort),wlists(iisource), &
          ifcharge,w(ichargesort), &
          ifdipole,w(idipstrsort),w(idipvecsort), &
          ifpot,w(ipot),iffld, w(ifld), &
          ntarget,w(itargetsort),wlists(iitarget), &
          ifpottarg,w(ipottarg),iffldtarg,w(ifldtarg), &
          epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp, &
          nboxes,laddr,nlev,scale,bsize,nterms, &
          wlists(iwlists),lwlists, nthread, &
	  ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
	  iffldsrc_far, iffldsrc_mpole, iffldsrc_near, &
	  ifpottarg_far, ifpottarg_mpole, ifpottarg_near, &
	  iffldtarg_far, iffldtarg_mpole, iffldtarg_near, &
	  potsrc_far, fldsrc_far, potsrc_mpole, fldsrc_mpole, &
	  potsrc_near, fldsrc_near, &
	  pottarg_far, fldtarg_far, pottarg_mpole, &
	  fldtarg_mpole, pottarg_near, fldtarg_near, ifnear_only)	

        t2=second()
!!$        t2=omp_get_wtime()
!        if( ifprint .eq. 1 ) call prin2('time in fmm main=*',t2-t1,1)

!       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0

        !if (ifprint .ge. 1) call prinf('lwlists=*',lused,1)
        !if (ifprint .ge. 1) call prinf('lused total =*',lused7,1)
       
        !if (ifprint .ge. 1) then
           !call prin2('memory / point = *',(lused7)/dble(nsource),1)
        !endif
	!call prin2('after w=*', w(1+lused7-100), 2*100)

!--------------------------------------------------------------------------
!	Reorder the arrays
!--------------------------------------------------------------------------

	! For the full FMM/ direct solver
        if(ifpot .eq. 1) then
          call l3dpsort(nsource,wlists(iisource),w(ipot),pot)
	endif

        if(iffld .eq. 1) then
	  !write(*,*) "yo"
          call l3dfsort(nsource,wlists(iisource),w(ifld),fld)
	endif

	! for the Multirate Solver
	!---------------------------------------------------
	! Namdi: For  multirate FMM
	! This is a problem because this is replacing the data 
	! that was calculated.
	! these functions use the info in wlists to write them in the 
	! proper order into the potential or field array
	!-------------------------------------------------

	if(ifnear_only .eq. 0) then

	  if(ifpottarg .eq. 1) then
	    call l3dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
	  endif
	  if(iffldtarg .eq. 1)  then
	    call l3dfsort(ntarget,wlists(iitarget),w(ifldtarg),fldtarg)
	  endif

	  ! multirate potential sources
	  if (ifpotsrc_far .eq. 1) then
	    call my_l3dpsort(nsource, wlists(iisource), potsrc_far, nthread)
	  endif

	  if (ifpotsrc_mpole .eq. 1) then
	    call my_l3dpsort(nsource, wlists(iisource), potsrc_mpole, nthread)
	  endif

	  if (ifpotsrc_near .eq. 1) then
	    call my_l3dpsort(nsource, wlists(iisource), potsrc_near, nthread)
	  endif

	  ! multirate field sources
	  if (iffldsrc_far .eq. 1) then
	    call my_l3dfsort(nsource, wlists(iisource), fldsrc_far, nthread)
	  endif

	  if (iffldsrc_mpole .eq. 1) then
	    call my_l3dfsort(nsource, wlists(iisource), fldsrc_mpole, nthread)
	  endif

	  if (iffldsrc_near .eq. 1) then
	    call my_l3dfsort(nsource, wlists(iisource), fldsrc_near, nthread)
	  endif

	  ! multirate potential targets
	  if (ifpottarg_far .eq. 1) then
	    call my_l3dpsort(ntarget,wlists(iitarget), pottarg_far, nthread)
	  endif

	  if (ifpottarg_mpole .eq. 1) then
	    call my_l3dpsort(ntarget,wlists(iitarget), pottarg_mpole, nthread)
	  endif

	  if (ifpottarg_near .eq. 1) then
	    call my_l3dpsort(ntarget,wlists(iitarget), pottarg_near, nthread)
	  endif

	  ! multirate field targets
	  if (iffldtarg_far .eq. 1) then
	    call my_l3dfsort(ntarget,wlists(iitarget), fldtarg_far, nthread)
	  endif

	  if (iffldtarg_mpole .eq. 1) then
	    call my_l3dfsort(ntarget,wlists(iitarget), fldtarg_mpole, nthread)
	  endif

	  if (iffldtarg_near .eq. 1) then
	    call my_l3dfsort(ntarget,wlists(iitarget), fldtarg_near, nthread)
	  endif

	endif !(ifnear_only .eq. 0)

	if (ifnear_only .eq. 1) then
	  
	  if(ifpotsrc_near .eq. 1) then
	    call l3dpsort(nsource,wlists(iisource),w(ipot),pot)
	    call my_l3dpsort(nsource, wlists(iisource), potsrc_near, nthread)
	  endif

	  if(iffldsrc_near .eq. 1) then
	    call my_l3dfsort(nsource, wlists(iisource), fldsrc_near, nthread)
	    !write(*,*) "Before sort, fld(1,1): ", fld(1,1)
	    call l3dfsort(nsource,wlists(iisource),w(ifld),fld)
	    !write(*,*) "After sort, fld(1,1): ", fld(1,1)
	  endif

	  if(ifpottarg_near .eq. 1) then
	    call l3dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
	    call my_l3dpsort(ntarget,wlists(iitarget), pottarg_near, nthread)
	  endif

	  if(iffldtarg_near .eq. 1)  then
	    call l3dfsort(ntarget,wlists(iitarget),w(ifldtarg),fldtarg)
	    call my_l3dfsort(ntarget,wlists(iitarget), fldtarg_near, nthread)
	  endif

	endif
	!----------------------------------------------------------	
	
        return
        

  end subroutine fmm

  !
  ! This is the direct solver
  !

  subroutine direct(nsource, source, ifcharge,charge, ifdipole, dipstr, dipvec, & 
			  ifpot, pot, iffld, fld, ntarget, target, ifpottarg, pottarg, &
			  iffldtarg, fldtarg, nthread)

	! Namdi: this is based on l3dpartdirect() of lfmm3dpart.f
	! Namdi: I added the nthread and made it able to handle an arbitrary amount
	! of threads
	! Namdi: Did some syntax editting so the code will work
	! Namdi: changed the declaring of function arguments
	! Namdi: added this
	use omp_lib

        implicit real *8 (a-h,o-z)

!       Laplace interactions in R^3: evaluate all pairwise particle
!       interactions (ignoring self-interaction) 
!       and interactions with targets via direct O(N^2) algorithm.
!
!       We use (1/r) for the Green's function,
!       without the (1/4 pi) scaling.  Self-interactions are not-included.
!   
!       INPUT PARAMETERS:
!
!       nsource: integer:  number of sources
!       source: real *8 (3,nsource):  source locations
!       ifcharge:  charge computation flag
!                  ifcharge = 1   =>  include charge contribution
!                                     otherwise do not
!       charge: complex *16 (nsource): charge strengths
!       ifdipole:  dipole computation flag
!                  ifdipole = 1   =>  include dipole contribution
!                                     otherwise do not
!       dipstr: complex *16 (nsource): dipole strengths
!       dipvec: real *8 (3,nsource): dipole orientation vectors. 
!
!       ifpot:  potential flag (1=compute potential, otherwise no)
!       iffld:  field flag (1=compute field, otherwise no)
!       ntarget: integer:  number of targets
!       target: real *8 (3,ntarget):  target locations
!       ifpottarg:  target potential flag 
!                   (1=compute potential, otherwise no)
!       iffldtarg:  target field flag 
!                   (1=compute field, otherwise no)
!
!       OUTPUT PARAMETERS:
!
!       pot: complex *16 (nsource): potential at source locations
!       fld: complex *16 (3,nsource): field (-gradient) at source locations
!       pottarg: complex *16 (ntarget): potential at target locations 
!       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
!


	!-----------------------------------------------
	! Namdi: I wrote the following code to pass in arguments
	!-------------------------------------------------
	
	! when the above line has "integer*4, intent(in), value :: ..." the code
	! starts missbehaving. For example, it will make nthread = -12344 and so on.
	integer, intent(in) :: nsource, ifcharge, ifdipole, ifpot, iffld, ifpottarg, &
					 iffldtarg, nthread

	real *8, intent(in) :: source(3, nsource), dipvec(3, nsource), target(3, ntarget)

	complex *16, intent(in) :: charge(nsource), dipstr(nsource)

	complex * 16, intent(out) :: pot(nsource), fld(3, nsource), pottarg(ntarget), &
				      fldtarg(3, ntarget)
	

        !dimension source(3,1),dipvec(3,1)
        !complex *16 charge(1),dipstr(1)
        !dimension target(3,1)
!
        !complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
!

	call omp_set_num_threads(nthread)

	! reset the potential and field 
	!if (ifpot .eq. 1) then
	  !call zero_scalar(pot, nsource, nthread)
	!endif
      
	!if (iffld .eq. 1) then
	  !call zero_vector(fld, nsource, nthread)
	!endif

	! zero out target potential and field
	!if (ifpottarg .eq. 1) then
	  !call zero_scalar(pottarg, ntarget, nthread)
	!endif

	!if (iffldtarg .eq. 1) then
	  !call zero_vector(fldtarg, ntarget, nthread)
	!endif

        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
        endif
        enddo
!       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        enddo
!
        if( ifpot .eq. 1 .or. iffld .eq. 1 ) then
	!$OMP PARALLEL DO DEFAULT(SHARED) &
	!$OMP PRIVATE(i,j,ptemp,ftemp) 
        do 6160 j=1,nsource
	  
	  ! Namdi: testing
	  !if (j .eq. 1) then
	    !write(*,*) "nthread: ", omp_get_num_threads()
	  !end if

        do 6150 i=1,nsource
            if (i .eq. j) goto 6150
            if (ifcharge .eq. 1 ) then
            call lpotfld3d(iffld,source(1,i),charge(i), &
                source(1,j),ptemp,ftemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (iffld .eq. 1) then
               fld(1,j)=fld(1,j)+ftemp(1)
               fld(2,j)=fld(2,j)+ftemp(2)
               fld(3,j)=fld(3,j)+ftemp(3)
            endif
            endif
            if (ifdipole .eq. 1) then
               call lpotfld3d_dp(iffld,source(1,i), &
                   dipstr(i),dipvec(1,i), &
                   source(1,j),ptemp,ftemp)
               if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
               if (iffld .eq. 1) then
                  fld(1,j)=fld(1,j)+ftemp(1)
                  fld(2,j)=fld(2,j)+ftemp(2)
                  fld(3,j)=fld(3,j)+ftemp(3)
               endif
            endif
 6150   continue
 6160   continue
	!$OMP END PARALLEL DO
        endif

        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
	!$OMP PARALLEL DO DEFAULT(SHARED) &
	!$OMP PRIVATE(i,j,ptemp,ftemp) 
        do j=1,ntarget
        do i=1,nsource
            if (ifcharge .eq. 1 ) then
            call lpotfld3d(iffldtarg,source(1,i),charge(i), &
                target(1,j),ptemp,ftemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (iffldtarg .eq. 1) then
               fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
               fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
               fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
            endif
            endif
            if (ifdipole .eq. 1) then
               call lpotfld3d_dp(iffldtarg,source(1,i), &
                   dipstr(i),dipvec(1,i), &
                   target(1,j),ptemp,ftemp) 
               if (ifpottarg .eq. 1 ) pottarg(j)=pottarg(j)+ptemp
               if (iffldtarg .eq. 1) then
                  fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                  fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                  fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
               endif
            endif
        enddo
        enddo
	!$OMP END PARALLEL DO
        endif
!
  end subroutine direct

  !
  ! This is the heart of the FMM
  !
  subroutine my_lfmm3dparttargmain( ier,iprec, ifevalfar,ifevalloc, &
    nsource, sourcesort, isource, ifcharge,chargesort, &
    ifdipole, dipstrsort, dipvecsort, ifpot, pot, iffld,& 
    fld, ntarget, targetsort, itarget, ifpottarg, pottarg, &
    iffldtarg, fldtarg, epsfmm, iaddr, rmlexp, mptemp, lmptemp, &
    nboxes, laddr, nlev, scale, bsize, nterms, wlists, lwlists, &
    nthread, ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near, &
    iffldsrc_far, iffldsrc_mpole, iffldsrc_near, ifpottarg_far, &
    ifpottarg_mpole, ifpottarg_near, iffldtarg_far, iffldtarg_mpole, iffldtarg_near,&
    potsrc_far, fldsrc_far, potsrc_mpole, &
    fldsrc_mpole, potsrc_near, fldsrc_near, &
    pottarg_far, fldtarg_far, pottarg_mpole, &
    fldtarg_mpole, pottarg_near,&
    fldtarg_near, ifnear_only )

    ! Namdi added this:
    use omp_lib

    implicit real *8 (a-h,o-z)

    integer, intent(in) 	:: isource(*), itarget(*), iaddr(2, nboxes), laddr(2, 200)
    real *8, intent(in)		:: sourcesort(3, *), targetsort(3,*), dipvecsort(3,*)
    complex *16, intent(in) 	:: chargesort(*), dipstrsort(*)        
    
    complex *16, intent(out) 	:: pot(*), fld(3,*), pottarg(*), fldtarg(3,*)    
    real *8, intent(out) 	:: wlists(*)
    
    real *8, intent(in) 	:: rmlexp(*)
    complex *16, intent(in) 	:: mptemp(lmptemp)
    dimension timeinfo(10)
    dimension center(3)    
    dimension scale(0:200)
    dimension bsize(0:200)
    dimension nterms(0:200)
    dimension list(10000)
    complex *16 		:: ptemp,ftemp(3)
    integer box(20)
    dimension center0(3),corners0(3,8)
    integer box1(20)
    dimension center1(3),corners1(3,8)
    dimension itable(-3:3,-3:3,-3:3)
    dimension wlege(40000)
    dimension nterms_eval(4,0:200)

    complex *16 :: ima
!	Namdi: added the following:
    integer :: nthread, ifnear_only
    integer :: ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near
    integer :: iffldsrc_far, iffldsrc_mpole, iffldsrc_near
    integer :: ifpottarg_far, ifpottarg_mpole, ifpottarg_near
    integer :: iffldtarg_far, iffldtarg_mpole, iffldtarg_near
	

    ! multirate FMM flags for the potential and fields
    integer :: ifpot_m, ifpottarg_m
    integer :: iffld_m, iffldtarg_m

    complex *16 :: potsrc_far(*), fldsrc_far(3,*) 
    complex *16 :: potsrc_mpole(*), fldsrc_mpole(3,*)
    complex *16 :: potsrc_near(*), fldsrc_near(3,*)

    complex *16 :: pottarg_far(*), fldtarg_far(3,*) 
    complex *16 :: pottarg_mpole(*), fldtarg_mpole(3,*)
    complex *16 :: pottarg_near(*), fldtarg_near(3,*)

    integer :: i, id, procs ! id is for testing

    data ima/(0.0d0,1.0d0)/

!
!
!     INPUT PARAMETERS:
!
!     iprec        precision flag (see above)   ELIMINATE???
!     ifevalfar    far field flag (1 means compute far field, 
!                                  else dont)
!     ifevalloc    local field flag (1 means compute local field, 
!                                    else dont)
!     nsource      number of sources
!     sourcesort   sorted source coordinates
!     isource      sorting index for sources
!     ifcharge     flag indicating potential includes contribution
!                  from charges
!     chargesort   sorted charge values
!     ifdipole     flag indicating potential includes contribution
!                  from dipoles
!     dipstrsort   sorted dipole strengths
!     dipvecsort   sorted dipole orientation vectors
!     ifpot        potential flag (1 => compute, else do not)
!     iffld        field flag (1 => compute, else do not)
!     ntarget      number of targets
!     targetsort   sorted array of target locations
!     itarget      sorting index for targets
!     ifpottarg    target potential flag (1 => compute, else do not)
!     iffldtarg    target field flag (1 => compute, else do not)
!     epsfmm       FMM tolerance
!     iaddr        iaddr(2,nboxes) array points to mpole/local
!                     expansions for each box
!     rmlexp       workspace to contain mpole/local expansions.
!     nboxes       number of boxes in FMM hierarchy
!     laddr        indexing array for FMM data structure
!     nlev         number of levels in FMM hierarchy
!     scale        array of scaling parameters
!     bsize        box dimension for FMM
!     nterms       array of nterms needed at each level
!     wlists       FMM data structure (real array)
!     lw           length of wlists
!
!
!     OUTPUT PARAMETERS:
!
!     pot          surface potential (if ifpot=1)
!     fld          surface field=-gradient(potential) (if iffld=1)
!     pottarg      target potential (if ifpot=1)
!     fldtarg      target field=-gradient(potential) (if iffld=1)
!     ier          error return code
!                  ier = 0    =>   normal execution
!                  ier = 4    =>   cannot allocate tree workspace
!                  ier = 8    =>   cannot alocate bulk FMM workspace
!                  ier = 16   =>   cannot allocate mpole exp workspace
!
!
!     ifprint is an internal information printing flag. 
!     Suppressed if ifprint=0.
!     Prints timing breakdown and other things if ifprint=1.
!     Prints timing breakdown, list information, and other things if ifprint=2.
!       
    !=========================================
    ! Namdi: Testing
    !========================================
    !write(*,*) "In lfmm3dparttargmain() in lfmm3dpart.f"

!   Namdi: I added this
    call omp_set_num_threads(nthread)

    ifprint=0
	
    ifpot_m=iunion(ifpotsrc_far,ifpotsrc_mpole,ifpotsrc_near)
    ifpottarg_m=iunion(ifpottarg_far,ifpottarg_mpole,ifpottarg_near)

    iffld_m=iunion(iffldsrc_far,iffldsrc_mpole,iffldsrc_near)
    iffldtarg_m=iunion(iffldtarg_far,iffldtarg_mpole,iffldtarg_near)

!	
!     
!    ... set the potential and field to zero
!
!    --------------------------
!    Namdi: Can make this parallel. Not optimized in Cache
!    ----------------------------
	
    !
    ! For the Full FMM
    !

    do i=1,nsource

      if( ifpot .eq. 1 .or. ifnear_only .eq. 1) then 
	pot(i)=0
      endif

      if( iffld .eq. 1 .or. ifnear_only .eq. 1) then
	fld(1,i)=0
	fld(2,i)=0
	fld(3,i)=0
      endif
    enddo

    do i=1,ntarget
      if( ifpottarg .eq. 1) pottarg(i)=0
      if( iffldtarg .eq. 1) then
	fldtarg(1,i)=0
	fldtarg(2,i)=0
	fldtarg(3,i)=0
      endif
    enddo

    !
    ! For the FMM near-field only
    !

    do i=1,nsource

      if(ifnear_only .eq. 1 .and. ifpotsrc_near .eq. 1) then 
	potsrc_near(i)=0
      endif

      if( ifnear_only .eq. 1 .and. iffldsrc_near .eq. 1) then
	fldsrc_near(1,i)=0
	fldsrc_near(2,i)=0
	fldsrc_near(3,i)=0
      endif

    enddo

    do i=1,ntarget

      if(ifnear_only .eq. 1 .and. ifpottarg_near .eq. 1) then 
	pottarg_near(i)=0
      endif

      if( ifnear_only .eq. 1 .and. iffldtarg_near .eq. 1) then
	fldtarg_near(1,i)=0
	fldtarg_near(2,i)=0
	fldtarg_near(3,i)=0
      endif

    enddo

    !
    ! For the FMM Separated 
    !
  
    do i=1,nsource

      if(ifpotsrc_near .eq. 1) potsrc_near(i)=0
      if(ifpotsrc_far .eq. 1) potsrc_far(i)=0

      if(iffldsrc_near .eq. 1) then
	fldsrc_near(1,i)=0
	fldsrc_near(2,i)=0
	fldsrc_near(3,i)=0
      endif

      if(iffldsrc_far .eq. 1) then
	fldsrc_far(1,i)=0
	fldsrc_far(2,i)=0
	fldsrc_far(3,i)=0
      endif

    enddo

    do i=1,ntarget

      if(ifpottarg_near .eq. 1) pottarg_near(i)=0
      if(ifpottarg_far .eq. 1) pottarg_far(i)=0

      if(iffldtarg_near .eq. 1) then
	fldtarg_near(1,i)=0
	fldtarg_near(2,i)=0
	fldtarg_near(3,i)=0
      endif

      if(iffldtarg_far .eq. 1) then
	fldtarg_far(1,i)=0
	fldtarg_far(2,i)=0
	fldtarg_far(3,i)=0
      endif

    enddo

    do i=1,10
      timeinfo(i)=0
    enddo

    !
    ! Start the far-field calculations
    !

    if( ifevalfar .eq. 0 ) goto 8000
       
!
!       ... initialize Legendre function evaluation routines
!
      nlege=100
      lw7=40000
      call ylgndrfwini(nlege,wlege,lw7,lused7)
!
      do i=0,nlev
	do itype=1,4
	  call l3dterms_eval(itype, epsfmm, nterms_eval(itype,i), ier)
        enddo
      enddo
!
      if (ifprint .ge. 2) then
	call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
      endif

      !=========================================
      ! Namdi: Testing
      !========================================      
!    
!       ... set all multipole and local expansions to zero
!
      !$omp parallel do default(shared) &
      !$omp private(ibox,box,center0,corner0,level, id, i_temp) &
      !$omp schedule(dynamic)
      do ibox = 1,nboxes
		
	call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))	
      enddo
      !$omp end parallel do         
      
      !=========================================
      ! Namdi: Testing
      !========================================

        if (ifprint .ge. 1) call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
!        t1=omp_get_wtime()
!
!       ... step 1, locate all charges, assign them to boxes, and
!       form multipole expansions
!
!!!        do 1200 ibox=1,nboxes
!!!        do 1300 ilev=3,nlev+1

      !$omp parallel do default(shared) &
      !$omp private(ibox,box,center0,corners0,level,npts,nkids,radius) &
      !$omp private(ier,i,j,ptemp,ftemp,cd) &
      !$omp schedule(dynamic)      

      do 1200 ibox=1,nboxes
!!!        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1	

	call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)

        level=box(1)

        if (ifprint .ge. 2) then
	  call prinf('ibox=*',ibox,1)
	  call prinf('box=*',box,20)
	  call prinf('nkids=*',nkids,1)
	endif

        if (nkids .eq. 0) then
	  ipts=box(14)
!	  npts=box(15)
!	  call prinf('ipts=*',ipts,1)
!	  call prinf('npts=*',npts,1)
	  npts=box(15)
	  if (ifprint .ge. 2) then
	    call prinf('npts=*',npts,1)
	    call prinf('isource=*',isource(box(14)),box(15))
	  endif
	endif
!
!       ... prune all sourceless boxes
!
	if( box(15) .eq. 0 ) goto 1200
!
        if (nkids .eq. 0) then
!
!       ... form multipole expansions
!
	  radius = (corners0(1,1) - center0(1))**2
	  radius = radius + (corners0(2,1) - center0(2))**2
	  radius = radius + (corners0(3,1) - center0(3))**2
	  radius = sqrt(radius)

	  call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))

	  if( ifcharge .eq. 1 ) then

            call l3dformmp_add_trunc(ier,scale(level), &
              sourcesort(1,box(14)),chargesort(box(14)), &
              npts,center0, &
              nterms(level),nterms_eval(1,level), &
              rmlexp(iaddr(1,ibox)),wlege,nlege)        
	  endif
 
	  if (ifdipole .eq. 1 ) then

            call l3dformmp_dp_add_trunc(ier,scale(level), &
              sourcesort(1,box(14)), &
              dipstrsort(box(14)),dipvecsort(1,box(14)), &
              npts,center0,nterms(level),nterms_eval(1,level), &
              rmlexp(iaddr(1,ibox)),wlege,nlege)            
	  endif
        endif
!
 1200    continue      
      !$omp end parallel do
      
 1300    continue

         t2=second()
!!$        t2=omp_get_wtime()
!!!        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
!       
        if (ifprint .ge. 1) call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
!!$        t1=omp_get_wtime()
!
!-----------------------------------------------------------------------
!       Step 2: In the adaptive FMM, a large leaf node may need to interact
!       with separated boxes at finer levels. This is called <list 3> in the
!       FMM. One takes individual sources in the large leaf node and 
!       maps them to local expansions in the target boxes.
!-----------------------------------------------------------------------
!      
      !$omp parallel do default(shared) &
      !$omp private(ibox,box,center0,corners0,level0,itype,list,nlist) &
      !$omp private(jbox,box1,center1,corners1,level1,ifdirect3,radius) &
      !$omp private(lused,ier,i,j,ptemp,ftemp,cd,ilist,npts)
      !!!$omp schedule(dynamic)

         do 3251 ibox=1,nboxes

         call d3tgetb(ier,ibox,box,center0,corners0,wlists)

!!!         if( box(15) .eq. 0 ) goto 3251

         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif

!       ... prune all sourceless boxes

!!!         if( box(15) .eq. 0 ) nlist=0


!       ... note that lists 3 and 4 are dual

!       ... form local expansions for all boxes in list 3
!       ... if target is childless, evaluate directly (if cheaper)
!====================================================================================
!	Maybe this should be calculated by the pottarg_far (farfield) no matter
!	if the calculation is cheaper or not.
!====================================================================================
       
!!!         call prinf('nlist3=*', nlist,1)
         do 3250 ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
        
            npts=box1(15)            
            if( npts .eq. 0 ) goto 3250

            level0=box(1)
            level1=box1(1)
!
!            ifdirect3 = 0
!            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
!     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
!
            ifdirect3 = 0

            if( ifdirect3 .eq. 0 ) then

               if( ifcharge .eq. 1 ) then

               call l3dformta_add_trunc(ier,scale(level0), &
                 sourcesort(1,box1(14)),chargesort(box1(14)), &
                 npts,center0, &
                 nterms(level0),nterms_eval(1,level0), &
                 rmlexp(iaddr(2,ibox)),wlege,nlege)

               endif

               if( ifdipole .eq. 1 ) then

               call l3dformta_dp_add_trunc(ier,scale(level0), &
                 sourcesort(1,box1(14)),dipstrsort(box1(14)), &
                 dipvecsort(1,box1(14)),npts,center0, &
                 nterms(level0),nterms_eval(1,level0), &
                 rmlexp(iaddr(2,ibox)),wlege,nlege)

               endif

            else

            call lfmm3dpart_direct_targ(box1,box,sourcesort, &
              ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
              ifpot,pot,iffld,fld, &
              targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)
	    !---------------------------------------------------------------
	    ! Namdi: Added this to store the far interactions 
	    !---------------------------------------------------------------
	
	    call lfmm3dpart_direct_targ(box1,box,sourcesort, &
              	ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
              	ifpotsrc_far,potsrc_far,iffldsrc_far, fldsrc_far, &
              	targetsort,ifpottarg_far,pottarg_far, &
     		iffldtarg_far,fldtarg_far)
	    endif
	    !----------------------------------------------------

            !endif
 3250    continue

 3251    continue
      !$omp end parallel do  
      

         t2=second()
!!$        t2=omp_get_wtime()
!!!        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
!
!-----------------------------------------------------------------------
!       Steps 3,4,5 are carried out by the routine lfmm3d_list2.
!       Step 3 is the merging of multipole expansions at every level.
!       Step 4 is the mapping of multipole expansions to local expansions
!              using <list 2>.
!       Step 5 is the recursive mapping of local expansions from parent to 
!              child.
!-----------------------------------------------------------------------
!
        if (ifprint .ge. 1) call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (iffld.eq.1) ifprune_list2 = 0

	if (ifpot_m .eq. 1) ifprune_list2 = 0
	if (iffld_m .eq. 1) ifprune_list2 = 0

        call lfmm3d_list2 (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm, &
	    timeinfo,wlists,mptemp,lmptemp,ifprune_list2)
!
!
!
        if (ifprint .ge. 1) call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
!!$        t1=omp_get_wtime()
!
!-----------------------------------------------------------------------
!       Step 6: In the adaptive FMM, a small leaf node may need to interact
!       with separated boxes at coarser levels. This is the dual of 
!       Step 2 and is called <list 4> in the FMM. 
!       The multipole expansion for the small leaf node is evaluated directly
!       at targets in the large target node.
!-----------------------------------------------------------------------
!
	!$omp parallel do default(shared) &
	!$omp private(ibox,box,center0,corners0,itype,list,nlist) &
	!$omp private(jbox,box1,center1,corners1,level1,ifdirect4,radius) &
	!$omp private(lused,ier,i,j,ptemp,ftemp,cd,ilist,level)
	!!!!!!$omp schedule(dynamic)
         do 3252 ibox=1,nboxes

         call d3tgetb(ier,ibox,box,center0,corners0,wlists)

!!!         if( box(15) .eq. 0 ) goto 3252

         itype=3
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif

!       ... prune all sourceless boxes

!!!         if( box(15) .eq. 0 ) nlist=0

!       ... note that lists 3 and 4 are dual


!       ... evaluate multipole expansions for all boxes in list 4 
!       ... if source is childless, evaluate directly (if cheaper)

!!!         call prinf('nlist4=*', nlist,1)
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)

            level=box1(1)

!            ifdirect4 = 0

!            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
!     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1

!           for future optimization - here, we just evaluate the 
!           multipole expansion, regardless of the number of sources
!           in the source box.

            ifdirect4 = 0

            if (ifdirect4 .eq. 0) then

!	    Namdi: Evaluate multipole expansion for the sources

	      if( box(15) .gt. 0 ) then
		call l3dmpevalall_trunc(scale(level),center1, &
              	rmlexp(iaddr(1,jbox)), &
              	nterms(level),nterms_eval(1,level), &
              	sourcesort(1,box(14)),box(15), &
              	ifpot,pot(box(14)), &
              	iffld,fld(1,box(14)), &
              	wlege,nlege,ier)
	      endif

!	    Namdi: Evaluate multipole expansion for the targets

	      if( box(17) .gt. 0 ) then
		call l3dmpevalall_trunc(scale(level),center1, &
              	rmlexp(iaddr(1,jbox)), &
              	nterms(level),nterms_eval(1,level), &
              	targetsort(1,box(16)),box(17), &
              	ifpottarg,pottarg(box(16)), &
              	iffldtarg,fldtarg(1,box(16)), &
              	wlege,nlege,ier)
	      endif

	    else !if( ifdirect4 .gt. 0)
            
	      call lfmm3dpart_direct_targ(box1,box,sourcesort, &
              ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
              ifpot,pot,iffld,fld, &
              targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)	    

            endif !THIS endif is for if(direct .eq. 4)

!---------------------------------------------------------------------
!234567	Namdi: I added this for testing the solver
!	changed the potentials and fields to the ones only designated
!	for the multipole interaction
!	lfmm3dpart_direct_targ changes the pot/field and pottarg/fldtarg
!	Redo the potential and field calculations with mypot_temp and myfld_temp
!	everytime the actual potential and field changes change the respective
!	temporary variables.
!------------------------------------------------------------------------------------
	    if (ifdirect4 .eq. 0) then

	      ! Namdi: Evaluate multirate multipole expansion for the sources
	      if( box(15) .gt. 0 ) then
		
               call l3dmpevalall_trunc(scale(level),center1, &
		rmlexp(iaddr(1,jbox)), &
		nterms(level),nterms_eval(1,level), &
		sourcesort(1,box(14)),box(15), &
		ifpotsrc_far, potsrc_far(box(14)), &
		iffldsrc_far, fldsrc_far(1,box(14)), &
		wlege,nlege,ier)
	      endif

	      ! Namdi: Evaluate multirate multipole expansion for the targets

		if( box(17) .gt. 0 ) then
              	 call l3dmpevalall_trunc(scale(level),center1, &
              	 rmlexp(iaddr(1,jbox)), &
              	 nterms(level),nterms_eval(1,level), &
              	 targetsort(1,box(16)),box(17), &
              	 ifpottarg_far,pottarg_far(box(16)), &
              	 iffldtarg_far,fldtarg_far(1,box(16)), &
              	 wlege,nlege,ier)
		endif

	    else ! if(iftest .eq. 1 .and. ifdirect4 .eq. 0)

	      ! Namdi: evaulate multirate direct calculation
	      
	      call lfmm3dpart_direct_targ(box1,box,sourcesort, &
              	 ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
              	 ifpotsrc_far, potsrc_far,iffldsrc_far, &
     		 fldsrc_far, targetsort,ifpottarg,pottarg_far, &
     		 iffldtarg, fldtarg_far)

	    endif !if (iftest .eq. 1. and. ifdirect4 .eq. 0)
!-------------------------------------------------------------------
        enddo
 3252   continue
	!$omp end parallel do	

        t2=second()
!!$        t2=omp_get_wtime()
!!!     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
!
!-----------------------------------------------------------------------
!       Step 7: Evaluate the local expansions for all relevant sources/targets.
!-----------------------------------------------------------------------
!
        if (ifprint .ge. 1) call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
!!$        t1=omp_get_wtime()

!       ... step 7, evaluate local expansions
!      and all fields directly

	!$omp parallel do default(shared) &
	!$omp private(ibox,box,center0,corners0,level,npts,nkids,ier)
	!!!!!!$omp schedule(dynamic)

        do 6201 ibox=1,nboxes

        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)

        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif

        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif

        if (nkids .eq. 0) then

!	------------------------------------------
!	Namdi: Evaluate local expansions (potential)
!	---------------------------------------------

!       ... evaluate local expansions
       
        level=box(1)
        npts=box(15)
       
        if (level .ge. 2) then

	  call l3dtaevalall_trunc(scale(level),center0, &
          rmlexp(iaddr(2,ibox)), &
          nterms(level),nterms_eval(1,level), &
          sourcesort(1,box(14)),box(15), &
          ifpot,pot(box(14)), &
          iffld,fld(1,box(14)), &
          wlege,nlege,ier)

	  call l3dtaevalall_trunc(scale(level),center0, &
          rmlexp(iaddr(2,ibox)), &
          nterms(level),nterms_eval(1,level), &
          targetsort(1,box(16)),box(17), &
          ifpottarg,pottarg(box(16)), &
          iffldtarg,fldtarg(1,box(16)), &
          wlege,nlege,ier)

! ------------------------------------------------
!	Namdi: for storing only the local expansion contribution	
!----------------------------------------------------
	    
	  call l3dtaevalall_trunc(scale(level),center0, &
          	rmlexp(iaddr(2,ibox)), &
          	nterms(level),nterms_eval(1,level), &
          	sourcesort(1,box(14)),box(15), &
          	ifpotsrc_far,potsrc_far(box(14)), &
          	iffldsrc_far,fldsrc_far(1,box(14)), &
          	wlege,nlege,ier)

	  call l3dtaevalall_trunc(scale(level),center0, &
         	rmlexp(iaddr(2,ibox)), &
          	nterms(level),nterms_eval(1,level), &
          	targetsort(1,box(16)),box(17), &
          	ifpottarg_far,pottarg_far(box(16)), &
          	iffldtarg_far,fldtarg_far(1,box(16)), &
          	wlege,nlege,ier)
!--------------------------------------------------------

!	for if (level. ge. 2) then
        endif
!	for if (nkids .eq. 9) then
        endif
!
 6201   continue
	!$omp end parallel do
	
        t2=second()
!!$        t2=omp_get_wtime()
!!!     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1


 8000   continue

	!
	! Start the near-field clculations
	!

        if( ifevalloc .eq. 0 ) goto 9000
 
        if (ifprint .ge. 1) call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
!!$        t1=omp_get_wtime()
!
!-----------------------------------------------------------------------
!       Step 8: Evaluate direct interactions locally
!-----------------------------------------------------------------------
!
	!$omp parallel do default(shared) &
	!$omp private(ibox,box,center0,corners0,nkids,list,nlist,npts) &
	!$omp private(jbox,box1,center1,corners1) &
	!$omp private(ier,ilist,itype)
	!!!!$omp schedule(dynamic)
        do 6202 ibox=1,nboxes

        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)

        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif

        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif

        if (nkids .eq. 0) then

!       ... evaluate self interactions

        if( box(15) .gt. 0 ) then
	   
	   ! calculate near-field interactions for regular FMM

     	   call lfmm3dpart_direct_self(box,sourcesort, &
          ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
          ifpot,pot,iffld,fld, &
          targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

!---------------------------------------------------------------
	!Namdi: To calculate the near-field interactions

	  call lfmm3dpart_direct_self(box,sourcesort, &
          	ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
          	ifpotsrc_near, potsrc_near,iffldsrc_near, &
     		fldsrc_near, targetsort,ifpottarg_near,pottarg_near, &
     	   	iffldtarg_near,fldtarg_near)
!---------------------------------------------------------------
	endif !if (box(15) .gt. 0)

!       ... retrieve list #1

!       ... evaluate interactions with the nearest neighbours

        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)

!       ... for all pairs in list #1, 
!       evaluate the potentials and fields directly

            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)

!       ... prune all sourceless boxes

         if( box1(15) .eq. 0 ) goto 6203


            call lfmm3dpart_direct_targ(box1,box,sourcesort, &
              ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
              ifpot,pot,iffld,fld, &
              targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

!---------------------------------------------------------------
!	Namdi: This is for evaluating the local interations that 
!	are solved directly
!	This function changes the field and potential of the 
!	sources.
!-------------------------------------------------------
	!Namdi: To calculate the near-field interactions

	    call lfmm3dpart_direct_targ(box1,box,sourcesort, &
              	ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
              	ifpotsrc_near, potsrc_near,iffldsrc_near, &
     		fldsrc_near, targetsort,ifpottarg_near, &
     		pottarg_near, iffldtarg_near,fldtarg_near)
!--------------------------------------------------------

 6203       continue
        endif

 6202   continue
	!$omp end parallel do
	
!!!        call prin2('inside fmm, pot=*',pot,2*nsource)


        t2=second()
!!$        t2=omp_get_wtime()
!!!     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1

 9000   continue

!!!        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)

        if (ifprint .ge. 1) call prin2('timeinfo=*',timeinfo,8)
       
        d=0
        do i=1,8
        d=d+timeinfo(i)
        enddo
       
        if (ifprint .ge. 1) call prin2('sum(timeinfo)=*',d,1)

        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('nsource=*',nsource,1)
        if (ifprint .ge. 1) call prinf('ntarget=*',ntarget,1)
       
        return
        end