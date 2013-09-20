!
! PFASST+OpenMM solution encapsulation.
!
! Each encapsulated LIBPFASST solution consists of two Fortran arrays:
!
!   * 'q' - position,
!   * 'p' - momentum.
!
! Namdi: used in Verlet scheme
!
module encap
  use iso_c_binding
  use pfasst
  implicit none

  type :: encap_t
     integer :: nparticles, ndim
     real(pfdp), pointer :: p(:, :), q(:, :), m(:)
  end type encap_t

  type :: encap_ctx_t
     integer :: nparticles, ndim
     real(pfdp), pointer :: m(:)
  end type encap_ctx_t

contains

  function get_m(ptr) result(r)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp), pointer :: r(:)

    type(encap_t), pointer :: sol
    call c_f_pointer(ptr, sol)

    r => sol%m
  end function get_m

  function get_p(ptr) result(r)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp), pointer :: r(:,:)

    type(encap_t), pointer :: sol
    call c_f_pointer(ptr, sol)

    r => sol%p
  end function get_p

  function get_q(ptr) result(r)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp), pointer :: r(:,:)

    type(encap_t), pointer :: sol
    call c_f_pointer(ptr, sol)

    r => sol%q
  end function get_q

  function nparticles(ptr) result(r)
    type(c_ptr), intent(in), value :: ptr
    integer :: r

    type(encap_t), pointer :: sol
    call c_f_pointer(ptr, sol)

    r = sol%nparticles
  end function nparticles

  function ndim(ptr) result(r)
    type(c_ptr), intent(in), value :: ptr
    integer :: r

    type(encap_t), pointer :: sol
    call c_f_pointer(ptr, sol)

    r = sol%ndim
  end function ndim

  subroutine encap_context_create(encap, nparticles, ndim, m)
    type(pf_encap_t), intent(out) :: encap
    integer,          intent(in)  :: nparticles, ndim
    real(pfdp),        intent(in) :: m(:)

    type(encap_ctx_t), pointer :: ctx
    allocate(ctx)

    ctx%nparticles = nparticles
    ctx%ndim       = ndim
    allocate(ctx%m(nparticles))
    ctx%m = m

    encap%ctx = c_loc(ctx)
    encap%create  => encap_create
    encap%destroy => encap_destroy
    encap%setval  => encap_setval
    encap%copy    => encap_copy
    encap%pack    => encap_pack
    encap%unpack  => encap_unpack
    encap%axpy    => encap_axpy
    encap%norm    => encap_norm
  end subroutine encap_context_create

  subroutine encap_context_destroy(encap)
    type(pf_encap_t),  intent(inout) :: encap
    type(encap_ctx_t), pointer :: ctx
    call c_f_pointer(encap%ctx, ctx)

    deallocate(ctx%m)
    deallocate(ctx)
  end subroutine encap_context_destroy

  subroutine encap_simple_create(sol, np, nd, m)
    type(c_ptr), intent(inout) :: sol
    integer,     intent(in)    :: np, nd
    real(pfdp),  intent(in)    ::  m(:)
    type(encap_t), pointer :: q

    allocate(q)

    q%nparticles = np
    q%ndim       = nd

    allocate(q%m(np))
    q%m = m

    allocate(q%q(np, nd))
    allocate(q%p(np, nd))

    sol = c_loc(q)
  end subroutine encap_simple_create

  ! Allocate/create solution (spatial data set) for the given level.
  !
  ! This is called for each SDC node.
  subroutine encap_create(sol, level, kind, nvars, shape, levelctx, encapctx)
    type(c_ptr),       intent(inout)     :: sol
    integer,           intent(in)        :: level, nvars, shape(:)
    integer,           intent(in)        :: kind
    type(c_ptr),       intent(in), value :: levelctx, encapctx

    type(encap_ctx_t), pointer :: ctx

    call c_f_pointer(encapctx, ctx)
    call encap_simple_create(sol, ctx%nparticles, ctx%ndim, ctx%m)
    

  end subroutine encap_create

  ! Deallocate/destroy solution.
  subroutine encap_destroy(ptr)
    type(c_ptr), intent(in), value :: ptr

    type(encap_t), pointer :: q
    call c_f_pointer(ptr, q)

    deallocate(q%m)
    deallocate(q%q)
    deallocate(q%p)
    deallocate(q)
  end subroutine encap_destroy

  ! Set solution value.
  subroutine encap_setval(ptr, val, flags)
    type(c_ptr), intent(in), value    :: ptr
    real(pfdp),  intent(in)           :: val
    integer,     intent(in), optional :: flags

    real(pfdp), pointer :: dst(:,:)
    integer :: which

    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (0)
       dst => get_q(ptr)
       dst = val
       dst => get_p(ptr)
       dst = val
    case (1)
       dst => get_p(ptr)
       dst = val
    case (2)
       dst => get_q(ptr)
       dst = val
    case default
       stop
    end select
  end subroutine encap_setval

  ! Copy solution value.
  subroutine encap_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in), value    :: dstptr, srcptr
    integer,     intent(in), optional :: flags

    real(pfdp), pointer :: dst(:,:), src(:,:)
    integer :: which

    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (0)
       dst => get_q(dstptr)
       src => get_q(srcptr)
       dst = src

       dst => get_p(dstptr)
       src => get_p(srcptr)
       dst = src
    case (1)
       dst => get_p(dstptr)
       src => get_p(srcptr)
       dst = src
    case (2)
       dst => get_q(dstptr)
       src => get_q(srcptr)
       dst = src
    case default
       stop
    end select

  end subroutine encap_copy

  ! Pack solution q into a flat array.
  subroutine encap_pack(z, ptr)
    type(c_ptr), intent(in), value  :: ptr
    real(pfdp),  intent(out)        :: z(:)

    integer :: n
    type(encap_t), pointer :: q

    call c_f_pointer(ptr, q)

    n = size(z) / 2
    z(:n)   = reshape(q%q, [ n ])
    z(n+1:) = reshape(q%p, [ n ])
  end subroutine encap_pack

  ! Unpack solution from a flat array.
  subroutine encap_unpack(ptr, z)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp),  intent(in)        :: z(:)
 
    integer :: n
    type(encap_t), pointer :: q

    call c_f_pointer(ptr, q)

    n = size(z) / 2
    q%q = reshape(z(:n),   [ q%nparticles, q%ndim ])
    q%p = reshape(z(n+1:), [ q%nparticles, q%ndim ])
  end subroutine encap_unpack

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine encap_axpy(yptr, a, xptr, flags)
    type(c_ptr), intent(in), value    :: xptr, yptr
    real(pfdp),  intent(in)           :: a
    integer,     intent(in), optional :: flags

    real(pfdp), pointer :: x(:,:), y(:,:)
    integer :: which

    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (0)
       x => get_p(xptr)
       y => get_p(yptr)
       y = a * x + y

       x => get_q(xptr)
       y => get_q(yptr)
       y = a * x + y
    case (1)
       x => get_p(xptr)
       y => get_p(yptr)
       y = a * x + y
    case (2)
       x => get_q(xptr)
       y => get_q(yptr)
       y = a * x + y
    case (12)
       x => get_p(xptr)
       y => get_q(yptr)
       y = a * x + y
    case default
       stop
    end select
  end subroutine encap_axpy

  function encap_norm(sol) result (norm)
    type(c_ptr), intent(in), value    :: sol
    real(pfdp) :: norm

    real(pfdp), pointer :: x(:, :), v(:, :)

    x => get_p(sol)
    v => get_q(sol)
    
    norm = max(maxval(abs(x)), maxval(abs(v)))
  end function encap_norm

  subroutine identity_interpolate(qF, qG, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qF, qG, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    real(pfdp), pointer :: dst(:, :), src(:, :)

    dst => get_q(qF)
    src => get_q(qG)
    dst = src

    dst => get_p(qF)
    src => get_p(qG)
    dst = src
  end subroutine identity_interpolate

  subroutine identity_restrict(qF, qG, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qF, qG, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    real(pfdp), pointer :: dst(:, :), src(:, :)

    dst => get_q(qG)
    src => get_q(qF)
    dst = src

    dst => get_p(qG)
    src => get_p(qF)
    dst = src
  end subroutine identity_restrict

end module encap
