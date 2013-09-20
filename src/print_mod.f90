module print_mod


contains

  subroutine print_cmplx(data, N)

    complex*16, intent(in) :: data(:,:)

    integer, intent(in) :: N

    integer :: i, j

    do i=1,N
      print *,  real( data(1:3, i) )
    enddo
  end subroutine

  subroutine print_real(data, N)

    real*8, intent(in) :: data(:,:)

    integer, intent(in) :: N

    integer :: i, j

    do i=1,N
      print *, data(1:3, i)
    enddo
  end subroutine

end module print_mod