c234567
	subroutine get_winfo(my_wlists, wlists, n)

	  implicit none

	  integer n, i

	  integer wlists(*), my_wlists(*)
 
	  !write(*,*) "In get_winfo() of my_func.f"
 
	  do i=1, n	    
	    my_wlists(i) = wlists(i)
	  enddo

	end subroutine
!------------------------------------------------------------------
	! Namdi: I created this function
	! Have to keep this function in the fortran 77 format
	! because it works. It doesn't for fortran 90.
	! I should parallelize this code
	subroutine my_l3dpsort(n,isource,pot, nthread)

	use omp_lib

        implicit real *8 (a-h,o-z)
        dimension isource(1)
        complex *16 pot(1)
	complex *16 temp(n)
	integer nthread
	
	call omp_set_num_threads(nthread)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
	
	! temp acts as psort
	do i=1, n
	temp(i) = pot(i)
	enddo
C$OMP END PARALLEL DO

c
ccc        call prinf('isource=*',isource,n)
c        
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
        do i=1,n
        pot(isource(i))=temp(i)
        enddo
C$OMP END PARALLEL DO
c
        return
        end	

	subroutine my_l3dfsort(n,isource,fld, nthread)

	use omp_lib

        implicit real *8 (a-h,o-z)

	! function arguments
        dimension isource(1)
        complex *16 fld(3,1)
	integer nthread
 
	! function variables
	complex *16 tempfld(3,n)

	call omp_set_num_threads(nthread)

ccc        call prinf('isource=*',isource,n)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
	do i=1, n
	  tempfld(1,i) = fld(1,i)
	  tempfld(2,i) = fld(2,i)
	  tempfld(3,i) = fld(3,i)
	enddo 
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
        do i=1,n
	  fld(1,isource(i))=tempfld(1,i)
	  fld(2,isource(i))=tempfld(2,i)
	  fld(3,isource(i))=tempfld(3,i)
        enddo
C$OMP END PARALLEL DO
c
        return
        end

c234567
	subroutine get_box_info(ier, ibox , box, center, corners, w)
	! Input: ibox, w
	! output: ier, box, center, conrners

	! this entry returns to the user the characteristics of
	! user-specified box  ibox.  

	  ! input parameters:

	  !  ibox - the box number for which the information is desired
	  !  w - storage area as created by the entry d3tstrcr (see above)

	  ! output parameters:

	  ! ier - the error return code.
	  ! ier=0 means successful execution
	  ! ier=2 means that ibox is either greater than the number of boxes
	  !           in the structure or less than 1.
	  ! box - an integer array dimensioned box(20). its elements describe 
	  !        the box number ibox, as follows:

	  ! 1. level - the level of subdivision on which this box 
	  ! 		was constructed; 
	  ! 2, 3, 4  - the coordinates of this box among  all
	  !   boxes on this level 
	      ! namdi: (x,y,z) coordinates? How is this an integer. 
	  ! Maybe 1 corresponds to 1 box over in a direction from an origin box?
	    ! The units are by boxes themselves

	  ! 5 - the daddy of this box, identified by it address
	  ! 		in array "boxes"
	  ! 6,7,8,9,10,11,12,13 - the  list of children of this box 
	      ! (eight of them, and the child is identified by its address
	      ! in the array "boxes"; if a box has only one child, only the
	      ! first of the four child entries is non-zero, etc.)

	  ! 14 - the location in the array iz of the particles 
		! living in this box
	  ! 15 - the number of particles living in this box
	  ! 16 - the location in the array iztarg of the targets
		! living in this box
	  ! 17 - the number of targets living in this box
	  ! 18 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
	  ! 19 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
	  ! 20 - reserved for future use

	  !  center - the center of the box number ibox 
	  !  corners - the corners of the box number ibox 

	  ! return to the user all information about the box ibox
c234567
	  implicit none

	  integer ibox, ier, box(20)

	  real *8 center(3), corners(3,8), w(*)

	  !write(*,*) "In get_box_info() of my_func.f"

	  call d3tgetb(ier, ibox, box, center, corners, w)
  
      end subroutine

c========================================================
c Namdi: I added the following function
c This function was added to avoid dealing with converting 
c things from .f files to .f90
c This function has support for an arbitrary amount of threads
c========================================================
	   subroutine my_lfmm3dparttargmain( ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists, nthread,
     $	   ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near,
     $	   iffldsrc_far, iffldsrc_mpole, iffldsrc_near,
     $	   ifpottarg_far, ifpottarg_mpole, ifpottarg_near,
     $	   iffldtarg_far, iffldtarg_mpole, iffldtarg_near,
     $	   potsrc_far, fldsrc_far, potsrc_mpole, 
     $	   fldsrc_mpole, potsrc_near, fldsrc_near,
     $	   pottarg_far, fldtarg_far, pottarg_mpole, 
     $	   fldtarg_mpole, pottarg_near,
     $	   fldtarg_near, iftest )

	! Namdi added this:
	use omp_lib

        implicit real *8 (a-h,o-z)
        dimension sourcesort(3,1), isource(1)
        complex *16 chargesort(1)
        complex *16 dipstrsort(1)
        dimension dipvecsort(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        dimension targetsort(3,1), itarget(1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        dimension wlists(1)
        dimension iaddr(2,nboxes)
        real *8 rmlexp(1)
        complex *16 mptemp(lmptemp)
        dimension timeinfo(10)
        dimension center(3)
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        dimension list(10 000)
        complex *16 ptemp,ftemp(3)
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
        dimension itable(-3:3,-3:3,-3:3)
        dimension wlege(40 000)
        dimension nterms_eval(4,0:200)
c
c	Namdi: added the following:
	integer nthread, iftest
	!iftarg_far, iftarg_mpole, iftarg_near
	integer ifpotsrc_far, ifpotsrc_mpole, ifpotsrc_near
	integer iffldsrc_far, iffldsrc_mpole, iffldsrc_near
	integer ifpottarg_far, ifpottarg_mpole, ifpottarg_near
	integer iffldtarg_far, iffldtarg_mpole, iffldtarg_near
	

	complex *16 potsrc_far(1), fldsrc_far(3,1) 
	complex *16 potsrc_mpole(1), fldsrc_mpole(3,1)
	complex *16 potsrc_near(1), fldsrc_near(3,1)

	complex *16 pottarg_far(1), fldtarg_far(3,1) 
	complex *16 pottarg_mpole(1), fldtarg_mpole(3,1)
	complex *16 pottarg_near(1), fldtarg_near(3,1)

	complex *16 myfld_temp(3,nsource)
	complex *16 mypot_temp(nsource)
	integer i

        data ima/(0.0d0,1.0d0)/

c
c
c     INPUT PARAMETERS:
c
c     iprec        precision flag (see above)   ELIMINATE???
c     ifevalfar    far field flag (1 means compute far field, 
c                                  else dont)
c     ifevalloc    local field flag (1 means compute local field, 
c                                    else dont)
c     nsource      number of sources
c     sourcesort   sorted source coordinates
c     isource      sorting index for sources
c     ifcharge     flag indicating potential includes contribution
c                  from charges
c     chargesort   sorted charge values
c     ifdipole     flag indicating potential includes contribution
c                  from dipoles
c     dipstrsort   sorted dipole strengths
c     dipvecsort   sorted dipole orientation vectors
c     ifpot        potential flag (1 => compute, else do not)
c     iffld        field flag (1 => compute, else do not)
c     ntarget      number of targets
c     targetsort   sorted array of target locations
c     itarget      sorting index for targets
c     ifpottarg    target potential flag (1 => compute, else do not)
c     iffldtarg    target field flag (1 => compute, else do not)
c     epsfmm       FMM tolerance
c     iaddr        iaddr(2,nboxes) array points to mpole/local
c                     expansions for each box
c     rmlexp       workspace to contain mpole/local expansions.
c     nboxes       number of boxes in FMM hierarchy
c     laddr        indexing array for FMM data structure
c     nlev         number of levels in FMM hierarchy
c     scale        array of scaling parameters
c     bsize        box dimension for FMM
c     nterms       array of nterms needed at each level
c     wlists       FMM data structure (real array)
c     lw           length of wlists
c
c
c     OUTPUT PARAMETERS:
c
c     pot          surface potential (if ifpot=1)
c     fld          surface field=-gradient(potential) (if iffld=1)
c     pottarg      target potential (if ifpot=1)
c     fldtarg      target field=-gradient(potential) (if iffld=1)
c     ier          error return code
c                  ier = 0    =>   normal execution
c                  ier = 4    =>   cannot allocate tree workspace
c                  ier = 8    =>   cannot alocate bulk FMM workspace
c                  ier = 16   =>   cannot allocate mpole exp workspace
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
	!=========================================
	! Namdi: Testing
	!========================================
	!write(*,*) "In lfmm3dparttargmain() in lfmm3dpart.f"

c	Namdi: I added this
	call omp_set_num_threads(nthread)

        ifprint=1

c	
c     
c       ... set the potential and field to zero
c
c	--------------------------
c	Namdi: Can make this parallel
c	----------------------------

	!Namdi: added setting mypot_temp and myfld_temp to zero
        do i=1,nsource

        if( ifpot .eq. 1) then 
	  pot(i)=0
	  if (iftest .eq. 1) then
	    mypot_temp(i) = 0
	  endif
	endif

        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
	   if (iftest .eq. 1) then
	    myfld_temp(1,i) = 0
	    myfld_temp(2,i) = 0
	    myfld_temp(3,i) = 0
	   endif
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        enddo

c
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
c
c       ... initialize Legendre function evaluation routines
c
        nlege=100
        lw7=40 000
        call ylgndrfwini(nlege,wlege,lw7,lused7)
c
        do i=0,nlev
        do itype=1,4
        call l3dterms_eval(itype,epsfmm,
     1       nterms_eval(itype,i),ier)
        enddo
        enddo
c
        if (ifprint .ge. 2) 
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
	
	!=========================================
	! Namdi: Testing
	!========================================
	!write(*,*) "set mpole and local expansions to zero"
c
c       ... set all multipole and local expansions to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,box,center0,corner0,level, id, i_temp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do ibox = 1,nboxes
	
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        enddo
C$OMP END PARALLEL DO
c
	!=========================================
	! Namdi: Testing
	!========================================
	!write(*,*) "Outside the do-loop"

        if (ifprint .ge. 1) call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
ccc        do 1200 ibox=1,nboxes
ccc        do 1300 ilev=3,nlev+1

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 1200 ibox=1,nboxes
ccc        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        level=box(1)
c
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
c        ipts=box(14)
c        npts=box(15)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(15)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(14)),box(15))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(15) .eq. 0 ) goto 1200
c
        if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
	    radius = (corners0(1,1) - center0(1))**2
	    radius = radius + (corners0(2,1) - center0(2))**2
	    radius = radius + (corners0(3,1) - center0(3))**2
	    radius = sqrt(radius)
c
            call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))

            if( ifcharge .eq. 1 ) then
c
            call l3dformmp_add_trunc(ier,scale(level),
     1         sourcesort(1,box(14)),chargesort(box(14)),
     $         npts,center0,
     $         nterms(level),nterms_eval(1,level),
     2         rmlexp(iaddr(1,ibox)),wlege,nlege)        
c
            endif
c 
            if (ifdipole .eq. 1 ) then

            call l3dformmp_dp_add_trunc(ier,scale(level),
     $         sourcesort(1,box(14)),
     1         dipstrsort(box(14)),dipvecsort(1,box(14)),
     $         npts,center0,nterms(level),nterms_eval(1,level),
     2         rmlexp(iaddr(1,ibox)),wlege,nlege)
            
            endif
         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c       
        if (ifprint .ge. 1) call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c-----------------------------------------------------------------------
c       Step 2: In the adaptive FMM, a large leaf node may need to interact
c       with separated boxes at finer levels. This is called <list 3> in the
c       FMM. One takes individual sources in the large leaf node and 
c       maps them to local expansions in the target boxes.
c-----------------------------------------------------------------------
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd,ilist,npts) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do 3251 ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
ccc         if( box(15) .eq. 0 ) goto 3251
c
         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
ccc         if( box(15) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
!====================================================================================
!	Maybe this should be calculated by the pottarg_far (farfield) no matter
!	if the calculation is cheaper or not.
!====================================================================================
c        
ccc         call prinf('nlist3=*', nlist,1)
         do 3250 ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            npts=box1(15)            
            if( npts .eq. 0 ) goto 3250
c
            level0=box(1)
            level1=box1(1)
c
c            ifdirect3 = 0
c            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
c     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
c
               if( ifcharge .eq. 1 ) then
c
               call l3dformta_add_trunc(ier,scale(level0),
     1            sourcesort(1,box1(14)),chargesort(box1(14)),
     $            npts,center0,
     $            nterms(level0),nterms_eval(1,level0),
     2            rmlexp(iaddr(2,ibox)),wlege,nlege)
c
               endif
c
               if( ifdipole .eq. 1 ) then

               call l3dformta_dp_add_trunc(ier,scale(level0),
     1            sourcesort(1,box1(14)),dipstrsort(box1(14)),
     2            dipvecsort(1,box1(14)),npts,center0,
     3            nterms(level0),nterms_eval(1,level0),
     $            rmlexp(iaddr(2,ibox)),wlege,nlege)

               endif
c
            else

            call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)
	    !---------------------------------------------------------------
	    ! Namdi: Added this to store the near interactions ( not farfield)
	    !------------------------------------------
	    if (iftest .eq. 1) then
	      ! test: mypot_temp, myfld_temp is being replaced
	      call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         	ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         	ifpotsrc_near,potsrc_near,iffldsrc_near, fldsrc_near,
     $         	targetsort,ifpottarg_near,pottarg_near,
     $		iffldtarg_near,fldtarg_near)
	    endif
	    !----------------------------------------------------

            endif
 3250    continue
c
 3251    continue
C$OMP END PARALLEL DO
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c-----------------------------------------------------------------------
c       Steps 3,4,5 are carried out by the routine lfmm3d_list2.
c       Step 3 is the merging of multipole expansions at every level.
c       Step 4 is the mapping of multipole expansions to local expansions
c              using <list 2>.
c       Step 5 is the recursive mapping of local expansions from parent to 
c              child.
c-----------------------------------------------------------------------
c
        if (ifprint .ge. 1) call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (iffld.eq.1) ifprune_list2 = 0
        call lfmm3d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,ifprune_list2)
c
c
c
        if (ifprint .ge. 1) call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c-----------------------------------------------------------------------
c       Step 6: In the adaptive FMM, a small leaf node may need to interact
c       with separated boxes at coarser levels. This is the dual of 
c       Step 2 and is called <list 4> in the FMM. 
c       The multipole expansion for the small leaf node is evaluated directly
c       at targets in the large target node.
c-----------------------------------------------------------------------
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd,ilist,level) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do 3252 ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
ccc         if( box(15) .eq. 0 ) goto 3252
c
         itype=3
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
ccc         if( box(15) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
ccc         call prinf('nlist4=*', nlist,1)
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box1(1)
c
c            ifdirect4 = 0
c
c            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
c     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
c           for future optimization - here, we just evaluate the 
c           multipole expansion, regardless of the number of sources
c           in the source box.
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then

c	    Namdi: Evaluate multipole expansion for the sources

	      if( box(15) .gt. 0 ) then
		call l3dmpevalall_trunc(scale(level),center1,
     $         	rmlexp(iaddr(1,jbox)),
     $         	nterms(level),nterms_eval(1,level),
     $         	sourcesort(1,box(14)),box(15),
     $         	ifpot,pot(box(14)),
     $         	iffld,fld(1,box(14)),
     $         	wlege,nlege,ier)
	      endif

c	    Namdi: Evaluate multipole expansion for the targets

	      if( box(17) .gt. 0 ) then
		call l3dmpevalall_trunc(scale(level),center1,
     $         	rmlexp(iaddr(1,jbox)),
     $         	nterms(level),nterms_eval(1,level),
     $         	targetsort(1,box(16)),box(17),
     $         	ifpottarg,pottarg(box(16)),
     $         	iffldtarg,fldtarg(1,box(16)),
     $         	wlege,nlege,ier)
	      endif

	    else !if( ifdirect4 .gt. 0)
            
	      call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)	    

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
	    if (iftest .eq. 1. .and. ifdirect4 .eq. 0) then

	      ! Namdi: Evaluate multirate multipole expansion for the sources
	      if( box(15) .gt. 0 ) then
		! Namdi: testing removing mypot_temp, myfld_temp
               call l3dmpevalall_trunc(scale(level),center1,
     $         rmlexp(iaddr(1,jbox)),
     $         nterms(level),nterms_eval(1,level),
     $         sourcesort(1,box(14)),box(15),
     $         ifpotsrc_mpole, potsrc_mpole(box(14)),
     $         iffldsrc_mpole, fldsrc_mpole(1,box(14)),
     $         wlege,nlege,ier)
	      endif

	      ! Namdi: Evaluate multirate multipole expansion for the targets

	      if (iftest .eq. 1) then

		if( box(17) .gt. 0 ) then
              	 call l3dmpevalall_trunc(scale(level),center1,
     $         	 rmlexp(iaddr(1,jbox)),
     $         	 nterms(level),nterms_eval(1,level),
     $         	 targetsort(1,box(16)),box(17),
     $         	 ifpottarg_mpole,pottarg_mpole(box(16)),
     $         	 iffldtarg_mpole,fldtarg_mpole(1,box(16)),
     $         	 wlege,nlege,ier)
		endif

	      endif ! if (iftest .eq. 1)
	    else ! if(iftest .eq. 1 .and. ifdirect4 .eq. 0)

	      ! Namdi: evaulate multirate direct calculation
	      ! myfld_temp, mypot_temp
	      call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         	 ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         	 ifpotsrc_near, potsrc_near,iffldsrc_near,
     $		 fldsrc_near, targetsort,ifpottarg,pottarg_near,
     $		 iffldtarg, fldtarg_near)

	    endif !if (iftest .eq. 1. and. ifdirect4 .eq. 0)
!-------------------------------------------------------------------
        enddo
 3252   continue
C$OMP END PARALLEL DO
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c
c-----------------------------------------------------------------------
c       Step 7: Evaluate the local expansions for all relevant sources/targets.
c-----------------------------------------------------------------------
c
        if (ifprint .ge. 1) call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6201 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
        if (nkids .eq. 0) then
c
c	------------------------------------------
c	Namdi: Evaluate local expansions (potential)
c	---------------------------------------------

c       ... evaluate local expansions
c       
        level=box(1)
        npts=box(15)
c       
        if (level .ge. 2) then

	  call l3dtaevalall_trunc(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),
     $     nterms(level),nterms_eval(1,level),
     $     sourcesort(1,box(14)),box(15),
     $     ifpot,pot(box(14)),
     $     iffld,fld(1,box(14)),
     $     wlege,nlege,ier)

	  call l3dtaevalall_trunc(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),
     $     nterms(level),nterms_eval(1,level),
     $     targetsort(1,box(16)),box(17),
     $     ifpottarg,pottarg(box(16)),
     $     iffldtarg,fldtarg(1,box(16)),
     $     wlege,nlege,ier)

c ------------------------------------------------
c	Namdi: for storing only the local expansion contribution
c	
c----------------------------------------------------
	  if (iftest .eq. 1 ) then
	    ! Namdi: testing with mypot_temp, myfld_temp
	    call l3dtaevalall_trunc(scale(level),center0,
     $     	rmlexp(iaddr(2,ibox)),
     $     	nterms(level),nterms_eval(1,level),
     $     	sourcesort(1,box(14)),box(15),
     $     	ifpotsrc_far,potsrc_far(box(14)),
     $     	iffldsrc_far,fldsrc_far(1,box(14)),
     $     	wlege,nlege,ier)

	    call l3dtaevalall_trunc(scale(level),center0,
     $    	rmlexp(iaddr(2,ibox)),
     $     	nterms(level),nterms_eval(1,level),
     $     	targetsort(1,box(16)),box(17),
     $     	ifpottarg_far,pottarg_far(box(16)),
     $     	iffldtarg_far,fldtarg_far(1,box(16)),
     $     	wlege,nlege,ier)
	  endif ! if (iftest .eq. 1)
c--------------------------------------------------------

c	for if (level. ge. 2) then
        endif
c	for if (nkids .eq. 9) then
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c 
        if (ifprint .ge. 1) call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c-----------------------------------------------------------------------
c       Step 8: Evaluate direct interactions locally
c-----------------------------------------------------------------------
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        if( box(15) .gt. 0 ) then

     	   call lfmm3dpart_direct_self(box,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,
     $     targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

!---------------------------------------------------------------
	!Namdi: To calculate the near-field interactions
	   if (iftest .eq. 1) then
		! Namdi: testing myfld_temp, mypot_temp
          	call lfmm3dpart_direct_self(box,sourcesort,
     $     	ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     	ifpotsrc_near, potsrc_near,iffldsrc_near,
     $		fldsrc_near, targetsort,ifpottarg_near,pottarg_near,
     $	   	iffldtarg_near,fldtarg_near)
	   endif
!---------------------------------------------------------------
	endif !if (box(15) .gt. 0)

c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and fields directly
c
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(15) .eq. 0 ) goto 6203
c
c
            call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

c---------------------------------------------------------------
c	Namdi: This is for evaluating the local interations that 
c	are solved directly
c	This function changes the field and potential of the 
c	sources.
c-------------------------------------------------------
	!Namdi: To calculate the near-field interactions

	    if (iftest .eq. 1) then
	      call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         	ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         	ifpotsrc_near, potsrc_near,iffldsrc_near,
     $		fldsrc_near, targetsort,ifpottarg_near,
     $		pottarg_near, iffldtarg_near,fldtarg_near)
	    endif
c--------------------------------------------------------
c
 6203       continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if (ifprint .ge. 1) call prin2('timeinfo=*',timeinfo,8)
c       
        d=0
        do i=1,8
        d=d+timeinfo(i)
        enddo
c       
        if (ifprint .ge. 1) call prin2('sum(timeinfo)=*',d,1)
c
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('nsource=*',nsource,1)
        if (ifprint .ge. 1) call prinf('ntarget=*',ntarget,1)
c       
        return
        end