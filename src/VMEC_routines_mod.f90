Module vmec_routines_mod
  Use kind_mod
  Use biotsavart ! from LIBSTELL
  Use coil_typedef
  Implicit None

  Public :: read_vmec_coils_file
  Public :: bfield_vmec_coils
  Public :: convert_vmec_coils_to_filaments

  Integer(int32), Parameter :: vmec_max_extcur = 100

  Integer(int32), Public :: vmec_nextcur = -1
  Real(real64), Public :: vmec_extcur(vmec_max_extcur) = 0.d0
  Type(coil_type), Public :: vmec_coil_new

  Private
  Save

Contains

  Subroutine read_vmec_coils_file(file_name,verbose)
    Implicit None
    Character(Len=*), Intent(In) :: file_name
    Logical, Intent(In) :: verbose

    If (verbose) Then
       Write(*,*) 'Reading vmec coils file:',Trim(Adjustl(file_name))
    End If

    Call parse_coils_file(Trim(file_name))
    vmec_nextcur= Size(coil_group)

    If (verbose) Then
       Write(*,*) 'vmec coils file had ',vmec_nextcur,' coil groups'
    End If

  End Subroutine read_vmec_coils_file

  Subroutine convert_vmec_coils_to_filaments
    Implicit None
    Integer(int32) :: ig, ic, npts_coil,npts_tot, i0, i1

    If (Sum(Abs(vmec_extcur)) .lt. 1.d-16) Then
       Write(*,*) 'Error, vmec_extcur not set!'
       Stop
    Endif

    npts_tot = 0
    Do ig = 1,vmec_nextcur      
       Do ic = 1,coil_group(ig)%ncoil
          npts_coil = Size(coil_group(ig)%coils(ic)%xnod,2)
          npts_tot = npts_tot + npts_coil
       Enddo
    Enddo

    Call allocate_coiltype_single(npts_tot,vmec_coil_new)
    vmec_coil_new%num_coils = 1
    vmec_coil_new%ncoil_pts = npts_tot
    i0 = 1
    i1 = 0
    Do ig = 1,vmec_nextcur
       Do ic = 1,coil_group(ig)%ncoil
          npts_coil = Size(coil_group(ig)%coils(ic)%xnod,2)
          i0 = 1 + i1
          i1 = i1 + npts_coil
          vmec_coil_new%coilxyz(i0:i1,1) = coil_group(ig)%coils(ic)%xnod(1,1:npts_coil)
          vmec_coil_new%coilxyz(i0:i1,2) = coil_group(ig)%coils(ic)%xnod(2,1:npts_coil)
          vmec_coil_new%coilxyz(i0:i1,3) = coil_group(ig)%coils(ic)%xnod(3,1:npts_coil)
          vmec_coil_new%current(i0:i1-1) = coil_group(ig)%coils(ic)%current*vmec_extcur(ig)
       Enddo
    Enddo

    !Write(*,*) 'Found npts_tot',npts_tot

  End Subroutine convert_vmec_coils_to_filaments


  Subroutine bfield_vmec_coils(r,phi,z,Npts,Bout,ierr)
    !   Bout = (:,[Br,Bz,Bt])
    Implicit None
    Integer(int32), Intent(In) :: Npts    
    Real(Real64), Intent(In), Dimension(Npts) :: r, z, phi
    Real(real64), Intent(Out), Dimension(Npts,3) :: Bout
    Integer(int32), Intent(Out) :: ierr
    ! Local variables
    Integer(int32) :: ig, numcoils, curindex, i
    Real(real64) :: br_tmp, bp_tmp, bz_tmp

    ierr = 0

    If (Sum(Abs(vmec_extcur)) .lt. 1.d-16) Then
       Write(*,*) 'Error, vmec_extcur not set!'
       Stop
    Endif

    Do i = 1,Npts
       Bout(i,1:3) = 0._real64

       Do ig = 1,vmec_nextcur
          numcoils = coil_group(ig)%ncoil
          curindex = Maxloc(Abs(coil_group(ig)%coils(1:numcoils)%current),1)
          br_tmp = 0.d0
          bp_tmp = 0.d0
          bz_tmp = 0.d0
          Call bfield(r(i),phi(i),z(i),br_tmp,bp_tmp,bz_tmp,ig)      
          Bout(i,1) = Bout(i,1) + br_tmp*vmec_extcur(ig)
          Bout(i,2) = Bout(i,2) + bz_tmp*vmec_extcur(ig)
          Bout(i,3) = Bout(i,3) + bp_tmp*vmec_extcur(ig)
       Enddo
    Enddo

  End Subroutine bfield_vmec_coils

End Module vmec_routines_mod
