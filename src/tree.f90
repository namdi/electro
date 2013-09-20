module tree

  use my_type

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! 
  ! Construct the tree object
  !
  
  subroutine tree_construct(obj, nsources, ntargets, iprec)
  
    integer, intent(in) :: nsources, ntargets, iprec

    type(tree_t), intent(out) :: obj
  
    integer :: i,j

    obj%nsources 	= nsources
    obj%ntargets 	= ntargets

    obj%iisource	= 0
    obj%iitarget	= 0
    obj%iwlists		= 0

    obj%nlev		= 0
    obj%nboxes		= 0
    obj%nbox		= 0
    obj%ntot		= 0

    obj%iprec		= iprec
    obj%epsfmm		= 0.0

    obj%lwlists		= 0
    obj%size		= 0
    obj%lused7		= 0

    ! zero out ladder
    do i=1,2
      do j=1,LADDER_MAX
	obj%ladder(i,j) = 0
      enddo
    enddo

    ! zero out center
    do i=1,3
      obj%center_temp(i) = 0.0
    enddo

    do j=1, 8
      do i=1,3
	obj%corners_temp(i,j) = 0.0
      enddo
    enddo

    nullify(obj%wlists)
    nullify(obj%boxes)
    nullify(obj%centers)
    nullify(obj%corners)

    allocate( obj%isource_tree(nsources) )
    allocate( obj%itarget_tree(ntargets) )

  end subroutine tree_construct

  !
  ! Destruct tree object
  !

  subroutine tree_destruct(obj)

    type(tree_t), intent(out) :: obj

    deallocate( obj%isource_tree )
    deallocate( obj%itarget_tree )

    deallocate( obj%wlists )
    deallocate( obj%boxes )
    deallocate( obj%centers )
    deallocate( obj%corners )

    nullify( obj%isource_tree )
    nullify( obj%itarget_tree )
    nullify( obj%wlists )
    nullify( obj%boxes )
    nullify( obj%centers )
    nullify( obj%corners )

  end subroutine tree_destruct

  !
  ! This function creates the tree
  !

  subroutine tree_create(tree, pos_src, pos_targ)

    real*8, intent(in) :: pos_src(:, :), pos_targ(:, :)

    type(tree_t), intent(out) :: tree

    real*8, pointer :: wlists(:)

    integer :: i, j, ier, nlev, nbox, ntot, iprec

    wlists	=> tree%wlists
    
    ! zero out ladder
    do i=1,2
      do j=1, LADDER_MAX
	tree%ladder(i,j) = 0
      enddo
    enddo

    ! zero out center_temp
    do i=1,3
     tree%center_temp(i) = 0.0
    enddo

    ! zero out corners_temp
    do j=1,8
      do i=1,3
	tree%corners_temp(i,j) = 0.0
      enddo
    enddo

    ! check if wilists is already allocated
    if ( associated(wlists) .eqv. .true.) then
      deallocate(wlists)
    endif
      
    ! get ntot (the length of wlists)
    call pre_create_tree(ier, tree%iprec, pos_src, tree%nsources, pos_targ, tree%ntargets, nlev, nbox, ntot)
    tree%ntot	= ntot

    ! allocate memory 
    allocate(tree%wlists(ntot))
    wlists => tree%wlists

    ! create the tree
    ier = 0
    call create_tree(ier, tree%iprec, pos_src, tree%nsources, pos_targ, tree%ntargets, tree%center_temp, tree%ladder, &
		      tree%nlev, tree%nboxes, tree%nbox, tree%epsfmm, tree%lused7, tree%iisource, tree%iwlists, &
		      tree%lwlists, tree%iitarget, tree%size, wlists, tree%ntot)

    ! set the boxes
    call set_boxes(tree)
    
  end subroutine tree_create

  !
  ! This routine is used to help create the FMM-tree.  It creates the box data.
  !

  subroutine set_boxes(tree)

  type(tree_t), intent(out) :: tree

  real*8, pointer	:: wlists(:)

  integer :: i, j, nboxes, nsources, ntargets, ier

  wlists	=> tree%wlists
  nboxes 	= tree%nboxes
  nsources	= tree%nsources
  ntargets	= tree%ntargets

  ! deallocate memory
  if( associated(tree%boxes) .eqv. .true. ) then
    deallocate(tree%boxes)
  endif
  
  if( associated(tree%centers) .eqv. .true. ) then
    deallocate(tree%centers)
  endif

  if( associated(tree%corners) .eqv. .true. ) then
    deallocate(tree%corners)
  endif

  ! assign memory to the following 
  allocate( tree%boxes(20, nboxes) )
  allocate( tree%centers(3, nboxes) )
  allocate( tree%corners(3, 8, nboxes) )

  ! fill iz(isource_tree) and iztarg(itarget_tree). These are the arrays that show where the target particles and the source
  ! particles are in the tree structure
  do i=1, nsources
    tree%isource_tree(i) = 0
  enddo

  do i=1, ntargets
    tree%itarget_tree(i) = 0
  enddo

  call get_winfo(tree%isource_tree, wlists(tree%iisource), nsources)
  call get_winfo(tree%itarget_tree, wlists(tree%iitarget), ntargets)

  ! zero out the following
  do i=1, 20
    tree%box(i) =0
  enddo

  do i=1,3
    tree%center_temp(i) = 0.0
  enddo

  do j=1,8
    do i=1,3
      tree%corners_temp(i,j) = 0.0
    enddo
  enddo

  ! store the information for all boxes into boxes
  ! i starts with index 1 in fortran
  ! i is the integer box number
  do i=1, nboxes
    call get_box_info(ier, i, tree%box, tree%center_temp, tree%corners_temp, wlists(tree%iwlists))
    call copy_box_data(tree, i)
  enddo

  end subroutine set_boxes

  !
  ! This routine is a helper routine that copies box data from 
  ! box to boxes, center_temp to centers, and corners_temp to corners
  !

  subroutine copy_box_data(tree, ibox)

    type(tree_t), intent(out) :: tree

    integer, intent(in) :: ibox

    integer :: i, j, k

    ! copy data from box to boxes(ibox)
    do i=1, 20
      tree%boxes(i, ibox) = tree%box(i)
    enddo

    ! copy data from center_temp to centers(ibox)
    do i=1,3
      tree%centers(i, ibox) = tree%center_temp(i)
    enddo

    ! copy data from corners_temp to corners(ibox)
    do j=1,8
      do i=1,3
	tree%corners(i,j,ibox) = tree%corners_temp(i,j)
      enddo
    enddo

  end subroutine copy_box_data

end module tree