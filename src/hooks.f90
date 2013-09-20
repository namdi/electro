!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none

 interface
     subroutine dump_mkdir(dname, dlen) bind(c)
       use iso_c_binding
       character(c_char), intent(in) :: dname
       integer(c_int),    intent(in), value :: dlen
     end subroutine dump_mkdir

     subroutine dump_numpy(dname, fname, endian, dim, shape, nvars, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in) :: dname, fname, endian(4)
       integer(c_int),    intent(in), value :: dim, nvars
       integer(c_int),    intent(in) :: shape(dim)
       real(c_double),    intent(in) :: array(nvars)
     end subroutine dump_numpy
  end interface

contains

  subroutine dump_hook(pf, level, state, ctx)
    use probin, only: outdir
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: ctx

    character(len=256)     :: fname
    type(ndarray), pointer :: qend

    call c_f_pointer(level%qend, qend)

    write(fname, "('s',i0.5,'i',i0.3,'l',i0.2,'.npy')") &
         state%step, state%iter, level%level

    call dump_numpy(trim(outdir)//c_null_char, trim(fname)//c_null_char, '<f8'//c_null_char, &
         qend%dim, qend%shape, size(qend%flatarray), qend%flatarray)

  end subroutine dump_hook

  !
  ! created my custum dump_hook to take in the output directory instead of using probin
  !
  subroutine my_dump_hook(pf, level, state, ctx)
    
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: ctx

    character(len=256)     :: fname
    type(ndarray), pointer :: qend
    
    call c_f_pointer(level%qend, qend)

    write(fname, "('s',i0.5,'i',i0.3,'l',i0.2,'.npy')") &
         state%step, state%iter, level%level

    call dump_numpy(trim(pf%ctx_char)//c_null_char, trim(fname)//c_null_char, '<f8'//c_null_char, &
         qend%dim, qend%shape, size(qend%flatarray), qend%flatarray)

  end subroutine my_dump_hook

!  subroutine echo_error(pf, level, state, ctx)
!    use iso_c_binding
!    use feval, only: exact
!    type(pf_pfasst_t), intent(inout) :: pf
!    type(pf_level_t),  intent(inout) :: level
!    type(pf_state_t),  intent(in)    :: state
!    type(c_ptr),       intent(in)    :: ctx

!    real(c_double) :: yexact(level%nvars)
!    real(pfdp), pointer :: qend(:)

!    qend => array1(level%qend)

!    call exact(state%t0+state%dt, yexact)
!    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
!         state%step+1, state%iter, maxval(abs(qend-yexact))
!  end subroutine echo_error


!  subroutine echo_residual(pf, level, state, ctx)
!    use iso_c_binding
!    use pf_mod_utils
!    type(pf_pfasst_t), intent(inout) :: pf
!    type(pf_level_t),  intent(inout) :: level
!    type(pf_state_t),  intent(in)    :: state
!    type(c_ptr),       intent(in)    :: ctx

!    real(pfdp), pointer :: r(:)

!    r => array1(level%R(level%nnodes-1))

!    print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
!         state%step+1, state%iter, level%level, maxval(abs(r))
!  end subroutine echo_residual


end module hooks
