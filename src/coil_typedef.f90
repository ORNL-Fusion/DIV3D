Module coil_typedef
  Use kind_mod, Only : real64, int32
  Implicit None

  Private

  ! If you add a component make sure to account for it in the
  ! allocate/deallocate routines below!
  Type, Public :: coil_type
    Real(real64), Allocatable :: coilxyz(:,:)  ! (ncoil_pts,3)
    Real(real64), Allocatable :: current(:)    ! (ncoil_pts)
    Integer(int32) :: ncoil_pts ! TOTAL number of points along created coil string
    Integer(int32) :: num_coils
  End Type coil_type

  Public :: allocate_coiltype
  Public :: allocate_coiltype_single
  Public :: deallocate_coiltype
  
Contains

  Subroutine allocate_coiltype_single(npts,coil)
    Implicit None
    Type(coil_type), Intent(Out) :: coil
    Integer(int32), Intent(In) :: npts
    
    Allocate(coil%coilxyz(npts,3))
    Allocate(coil%current(npts))
    coil%ncoil_pts = npts
    coil%num_coils = 1
    coil%coilxyz = 0._real64
    coil%current = 0._real64
  End Subroutine allocate_coiltype_single    
  
  Subroutine allocate_coiltype(num_coils,ntorpts,coil)
    Implicit None
    Type(coil_type), Intent(Out) :: coil
    Integer(int32), Intent(In) :: num_coils, ntorpts
    Integer(int32) :: npts
    npts = num_coils*(ntorpts*2 + 1)   
    Allocate(coil%coilxyz(npts,3))
    Allocate(coil%current(npts))
    coil%ncoil_pts = npts
    coil%num_coils = num_coils
    coil%coilxyz = 0._real64
    coil%current = 0._real64
  End Subroutine allocate_coiltype

  Subroutine deallocate_coiltype(coil)
    Implicit None
    Type(coil_type), Intent(inout) :: coil
    Deallocate(coil%coilxyz,coil%current)
    coil%ncoil_pts = 0
    coil%num_coils = 0
  End Subroutine deallocate_coiltype
  
End Module coil_typedef
