module direct
  use iso_c_binding		! makes the fortran code look like c
contains

  subroutine calculate_field(field, potential, source, nsources, charge, fieldtarg, pottarg, &
			      target, ntargets, nthreads) bind(c, name="calculate_field")
    
    use omp_lib
    
    implicit none

    ! Function Arguments
    integer(c_int), intent(in), value :: nsources, nthreads, ntargets

    real(c_double), intent(in) :: source(3, nsources), charge(nsources), target(3, ntargets)

    real(c_double), intent(out) :: field(3, nsources), potential(nsources), fieldtarg(3, ntargets), &
				   pottarg(ntargets)

    ! Local variables
    real(c_double) :: dist, distSq, dx, dy, dz, eMag, ex, ey, ez, qi, qj, pot, K

    integer(c_int) :: i, j

    call omp_set_num_threads(nthreads)

    K = 1.0

    !$omp parallel do &
    !$omp shared(K, source, charge, nsources, potential, field) &
    !$omp private(dist, distSq, dx, dy, dz, eMag, ex, ey, ez, qi, qj, pot, i, j) &
    !$omp& schedule(dynamic)

    do i = 1, nsources
      ! reset the E-field variables
      eMag = 0.0
      ex = 0.0
      ey = 0.0
      ez = 0.0
      pot = 0.0

      qi = charge(i)

      do j = 1, nsources
	if (i .ne. j) then
	    ! set the charge for the jth particle
	    qj = charge(j)

	    ! calculate the distance
	    dx = source(1,j) - source(1,i)
	    dy = source(2,j) - source(2,i)
	    dz = source(3,j) - source(3,i)
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
      field(1,i) = ex
      field(2,i) = ey
      field(3,i) = ez

      ! store the potential for the ith particle
      potential(i) = pot
    end do
    !$omp end parallel do

    !$omp parallel do &
    !$omp shared(K, source, target, charge, nsources, ntargets, pottarg, fieldtarg) &
    !$omp private(dist, distSq, dx, dy, dz, eMag, ex, ey, ez, qi, qj, pot, i, j) &
    !$omp& schedule(dynamic)
    do i = 1, ntargets
      pot  = 0.0
      eMag = 0.0
      ex   = 0.0
      ey   = 0.0
      ez   = 0.0
      qi   = 1.0

      do j = 1, nsources
	! the charge for the jth source
	qj = charge(j)

	! calculate the distance
	dx = source(1,j) - target(1,i)
	dy = source(2,j) - target(2,i)
	dz = source(3,j) - target(3,i)
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
      fieldtarg(1,i) = ex
      fieldtarg(2,i) = ey
      fieldtarg(3,i) = ez

      ! store the potential at the ith targtet
      pottarg(i) = pot
     
    enddo
    !$omp end parallel do

  end subroutine calculate_field

  subroutine forwardEuler(pos, vel, nparticles, dt, nthreads) bind(c, name = "forwardEuler")
    
    use omp_lib

    implicit none

    integer(c_int), intent(in), value :: nparticles, nthreads

    real(c_double), intent(in), value :: dt

    real(c_double), intent(in) :: vel(3, nparticles)

    real(c_double), intent(out) :: pos(3, nparticles)

    !integer(c_int) :: i

    integer :: i

    write(*,*) "In forwardEuler() of direct.f90"
    
    call omp_set_num_threads(nthreads)

    !$omp parallel do &
    !$omp shared(nparticles, dt, pos, vel) &
    !$omp private(i) &
    !$omp schedule(dynamic)

    do i = 1, nparticles
      pos(1,i) = pos(1,i) + vel(1,i)*dt
    end do
    !$omp end parallel do

    !$omp parallel do &
    !$omp shared(nparticles, dt, pos, vel) &
    !$omp private(i) &
    !$omp schedule(dynamic)

    do i = 1, nparticles
      pos(2,i) = pos(2,i) + vel(2,i)*dt
    end do
    !$omp end parallel do

    !$omp parallel do &
    !$omp shared(nparticles, dt, pos, vel) &
    !$omp private(i) &
    !$omp schedule(dynamic)

    do i = 1, nparticles
      pos(3,i) = pos(3,i) + vel(3,i)*dt
    end do
    !$omp end parallel do

  end subroutine forwardEuler
end module direct