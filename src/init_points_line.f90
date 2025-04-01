Module initialize_points
  Implicit None
  Public :: init_points_line

Contains
  !-----------------------------------------------------------------------------
  !+ 
  !-----------------------------------------------------------------------------
  Subroutine init_points_line
    ! Inputs: 
    ! Author(s): J.D. Lore - 07/26/2011 - xxx

    ! Modules used:
    Use kind_mod, Only : int32, real64
    Use io_unit_spec, Only: iu_launch, iu_surf
    Use run_settings_namelist, Only : period, fname_surf,npts_start,fname_launch
    Use math_routines_mod, Only : wrap_phi

    Implicit none

    Integer(int32) :: ii, rand_ind
    Integer(int32) :: npts_line,nip0, ip_step
    Real(real64) :: P1, rnum
    Real(real64),Allocatable,Dimension(:) :: rsurf,zsurf,phisurf

    !- End of header -------------------------------------------------------------

    ! Read surface file
    Write(*,*) 'Reading surface data from ',Trim(Adjustl(fname_surf))
    Open(iu_surf,file=fname_surf)
    Read(iu_surf,*) period, nip0, ip_step
    Read(iu_surf,*) npts_line
    Allocate(rsurf(npts_line),zsurf(npts_line),phisurf(npts_line))
    Do ii = 1,npts_line
       Read(iu_surf,*) rsurf(ii)
       Read(iu_surf,*) zsurf(ii)
       Read(iu_surf,*) phisurf(ii)
    Enddo
    Close(iu_surf)

    Write(*,*) 'Writing launch point data to ',Trim(Adjustl(fname_launch))
    Open(iu_launch,file=fname_launch)
    Write(iu_launch,*) npts_start

    Do ii = 1,npts_start
       ! Choose an integer between 1 and npts_line
       Call Random_number(rnum)
       rand_ind = Nint(npts_line*rnum)
       P1 = phisurf(rand_ind)
       Call wrap_phi(P1,period)
       Write(iu_launch,*) rsurf(rand_ind),zsurf(rand_ind),P1
    Enddo
    Close(iu_launch)

  End Subroutine init_points_line
  !-----------------------------------------------------------------------------

End Module initialize_points
