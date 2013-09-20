module probin
  use pf_mod_dtype

  character(len=64), save :: outdir ! directory name for output
  

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename

    integer :: un

    namelist /prbin/ &
         outdir


    !
    ! defaults
    !

    outdir       = ""

    !
    ! read
    !

    un = 66
    open(unit=un, file=filename, status='old', action='read')
    read(unit=un, nml=prbin)
    close(unit=un)


  end subroutine probin_init

end module probin

