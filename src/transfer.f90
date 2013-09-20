!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer spatial (interpolate, restrict) routines.

module transfer
  use feval
  
  ! testing
  use print_mod
  implicit none
contains

  !
  ! interpolate the spaital data
  !

  subroutine interpolate(qFp, qGp, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    real(pfdp),      pointer :: qF3(:,:,:), qG3(:,:,:)

    integer :: N
    
    ! qF3 and qG3 are of size [3, N, 2]
    qF3	=> array3(qFp)
    qG3	=> array3(qGp)

    N = size(qF3, 2)

    ! idenity interpolator (simply copies) because this is a particle/ gridless method
    qF3(1:3, 1:N, 1:2) = qG3(1:3, 1:N, 1:2)

  end subroutine interpolate

  !
  ! restrict spatially
  !

  subroutine restrict(qFp, qGp, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    real(pfdp), pointer :: qF3(:,:,:), qG3(:,:,:)

    integer :: N

    ! qF3 is of size [3, N, 2]
    qF3	=> array3(qFp)
    qG3	=> array3(qGp)

    N = size(qF3, 2)

    ! idenity restrictor (simply copies) because this is a particle/ gridless method
    qG3(1:3, 1:N, 1:2) = qF3(1:3, 1:N, 1:2)

  end subroutine restrict

end module transfer
