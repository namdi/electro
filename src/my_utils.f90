!
! how to read parameters from a file
!
module my_utils
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none

contains

  !===============================================================
  !===============================================================

  !
  ! load parameters
  !
  subroutine load_parameters(fname, nsources, ntargets, distribution, dt, nsteps)


    character(len=*), intent(in) :: fname

    integer, intent(out) :: nsources, ntargets, distribution, nsteps

    double precision, intent(out) :: dt

    character(len=300) :: buffer, temp

    integer :: f1, idx, length  

    length = len_trim(fname)
    f1 = 100

    open(f1, file=trim(fname), status="old")

    ! store nsources
    read(f1, '(A)') buffer
    idx = len("nsources: ")
    length = len_trim(buffer)
    read( buffer(idx:length) ,* ) nsources
  
    ! store ntargets
    read(f1, '(A)') buffer
    idx = len("ntargets: ")
    length = len_trim(buffer)
    read( buffer(idx:length), * ) ntargets

    ! store distribution
    read(f1, '(A)') buffer
    idx = len("distribution: ")
    length = len_trim(buffer)
    read( buffer(idx:length), * ) distribution
  
    ! store dt
    read(f1, '(A)') buffer
    idx = len("dt: ")
    length = len_trim(buffer)
    read( buffer(idx:length), * ) dt

    ! store nsteps: PROBLEM IS HERE
    read(f1, '(A)') buffer
    idx = len("nsteps: ")
    length = len_trim(buffer)
    read( buffer(idx:length), * ) nsteps  

    ! close file
    close(f1)

  end subroutine load_parameters


  ! 
  ! load the charge 
  !

  subroutine mu_load_charge(fname, charge, N)

    integer, intent(in)		:: N
  
    complex*16, intent(out) 	:: charge(N)

    integer :: i, u

    character(len=300) :: fname

    real*8 :: temp(2, N)

      u = 100

      open(u, file=fname, status='old')

      do i=1, N
	read(u, *, end=110) temp(1,i), temp(2,i)
	charge(i) = cmplx( temp(1,i), temp(2,i) )
      enddo

110   close(u)

  end subroutine mu_load_charge

  !
  ! load dipole information
  !

  subroutine mu_load_dipole(fname_str, fname_vec, dipstr, dipvec, N)
    
    character(len=*), intent(in):: fname_str, fname_vec

    integer, intent(in)		:: N
  
    complex*16, intent(out) 	:: dipstr(N)

    real*8, intent(out)		:: dipvec(3,N)

    integer 	:: i, u

    real*8	:: temp(2,N)

      u = 100
   
      open(u, file=fname_str, status='old')

      do i=1, N
	read(u, *, end=210) temp(1, i), temp(2, i)
	dipstr(i) = cmplx(temp(1,i), temp(2,i))
      enddo

210   close(u)

      open(u, file=fname_vec, status='old')

      do i=1, N
	read(u, *, end=220) dipvec(1,i), dipvec(2,i), dipvec(3,i)
      enddo

220   close(u)

  end subroutine mu_load_dipole

  !
  !
  ! load the initial condition
  !

  subroutine mu_load_ic(fdir, q, nsources, ntargets)
    

    character(len=*), intent(in) :: fdir

    integer, intent(in) :: nsources, ntargets

    real*8, intent(out) :: q(3, nsources + ntargets, 2)

    integer :: f1, N, length, j
   
    character(len=300) :: fname, buffer
    
    N = nsources + ntargets
    
    ! read the source positions

    f1 = 100
    length = len_trim(fdir)     
    fname = fdir(1:length) // "/pos_source.txt"    
    
    open(f1, file=fname, status='old')

    do j=1,nsources
      read(f1, *, end=10) q(1,j, 1), q(2,j, 1), q(3,j, 1)  
    enddo

10  close(f1)

    ! read the source velocities

    fname = fdir(1:length) // "/vel_source.txt"

    open(f1, file=fname, status='old')    

    do j=1,nsources
      read(f1, *, end=20) q(1,j,2), q(2,j,2), q(3,j,2)
    enddo

20    close(f1)
    
    ! read the target positions

    fname = fdir(1:length) // "/pos_target.txt"

    open(f1, file=fname, status='old')    

    do j=nsources+1,N
      read(f1, *, end=30) q(1,j,1) , q(2,j,1), q(3,j,1)
    enddo

30    close(f1)

    ! read the target velocities

    fname = fdir(1:length) // "/vel_target.txt"

    open(f1, file=fname, status='old')    

    do j=nsources+1,N
      read(f1, *, end=40) q(1,j,2), q(2,j,2), q(3,j,2)
    enddo

40    close(f1)
    
  end subroutine mu_load_ic

end module my_utils