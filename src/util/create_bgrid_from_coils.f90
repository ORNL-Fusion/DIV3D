Program create_bgrid_from_coils
  !
  ! Author(s): J.D. Lore
  !
  Use kind_mod, Only : int32, real64
  Use create_bgrid_from_coils_namelist
  Use initialize_bfield_div3d, Only : init_bfield
  Use math_routines_mod, Only : rlinspace
  Use setup_bfield_module, Only : nfp_bfield, bfield
  Use bfield_module, Only : calc_B_rzphi_general
  Use phys_const, Only : pi
  Use parallel_mod
  Use timing_mod, Only : get_elapsed_time, init_timing  
  Implicit None
  
  ! Local variables
  Logical :: verbose
  Integer :: tstart ! Timing variable
  Real(real64) :: dr, dz, dphi
  Real(real64), Allocatable :: R(:), Z(:), P(:)
  Real(real64), Allocatable :: Br(:,:,:), Bz(:,:,:), Bphi(:,:,:)
  Integer(int32) :: ierr_out, ir, jz, kt
  

  !---------------------------------------------------------------------------
  ! 0. Initialization
  !---------------------------------------------------------------------------
  Call init_mpi
  Call init_timing

  verbose = .false.
  If (rank == 0) Then
     Call system_clock(tstart)
     verbose = .true.
     Write(*,'(/A)') '-------------------------------------------------------------------------'
     Write(*,'(a)') " Starting create_bgrid_from_coils"
  End If

  ! Read namelists
  Call read_create_bgrid_from_coils_namelist(verbose)
  
  !----------------------------------------------------------
  ! 1. Initialize magnetic field
  !----------------------------------------------------------  
  Call init_bfield(verbose)

  !----------------------------------------------------------
  ! 2. Allocate variables and setup grid
  !----------------------------------------------------------
    
  dr = (Rmax-Rmin)/(nr-1)
  dz = (Zmax-Zmin)/(nz-1)
  dphi = 2.d0*pi/nfp_bfield/(nphi-1)

  Write(*,*) "Creating bgrid using the following settings:"

  Write(*,*) "Prefix: ",Trim(outfile_prefix)
  
  Write(*,*) '[Rmin,Rmax] = ',Rmin,Rmax,' m'
  Write(*,*) '[Zmin,Zmax] = ',Zmin,Zmax,' m'
  Write(*,*) '[Phimin,Phimax] = ',0,period*180.d0/pi,' deg'

  Write(*,*) '[nR,nZ,nPhi] = ',nr,nz,nphi
  Write(*,*) 'Note: nPhi is number of planes including 0 and 2*pi/nfp_bfield'
  Write(*,*) 'The last (redundant) plane is not written to the grid file'
  
  Write(*,*) 'dR   = ',dr,' (m)'
  Write(*,*) 'dZ   = ',dz,' (m)'
  Write(*,*) 'dphi = ',dphi*180.d0/pi,' (deg)'

  Allocate(R(nr))
  Allocate(Z(nz))
  Allocate(P(nphi))
  
  R = rlinspace(Rmin,Rmax,nr)
  Z = rlinspace(Zmin,Zmax,nz)
  P = rlinspace(0.d0,period,nphi)

  Allocate(Br(nr,nz,nphi-1))
  Allocate(Bz(nr,nz,nphi-1))
  Allocate(Bphi(nr,nz,nphi-1))

  !----------------------------------------------------------
  ! 3. Evaluate B
  !----------------------------------------------------------
  
  Do kt = 1,nphi
     If (verbose) Write(*,*) 'Working on slice ',kt,' of ',nphi
     Do jz = 1,nz
        Do ir = 1,nr
           Call calc_B_rzphi_general(bfield,R(ir),Z(jz),P(kt),1, &
                br(ir,jz,kt),bz(ir,jz,kt),bphi(ir,jz,kt),ierr_out)

           If (ierr_out .ne. 0) Then
              Write(*,*) 'Bfield error at [ir,jz,kt] = ',ir,jz,kt
              Write(*,*) 'Error flag:',ierr_out
              Call fin_mpi(.true.)
           End If
        End Do
     End Do
  End Do

  !----------------------------------------------------------
  ! 4. Write grid
  !----------------------------------------------------------
  Open(unit=99, file=trim(outfile_prefix) // '_layout.dat', status='replace', action='write')
  Write(99, '(4I4,4F10.5)') nr, nz, nphi-1, nfp_bfield, Rmin, Rmax, Zmin, Zmax
  Close(99)

  Write(*,*) 'Writing Br'
  Open(unit=99, file=trim(outfile_prefix) // '_r.dat', status='replace', action='write')
  Do ir = 1,nr
     Do jz = 1,nz
        Do kt = 1,nphi-1
           Write(99,'(E18.12)') Br(ir,jz,kt)
        End Do
     End Do
  End Do
  Close(99)

  Write(*,*) 'Writing Bz'
  Open(unit=99, file=trim(outfile_prefix) // '_z.dat', status='replace', action='write')
  Do ir = 1,nr
     Do jz = 1,nz
        Do kt = 1,nphi-1
           Write(99,'(E18.12)') Bz(ir,jz,kt)
        End Do
     End Do
  End Do
  Close(99)

  Write(*,*) 'Writing Bphi'
  Open(unit=99, file=trim(outfile_prefix) // '_phi.dat', status='replace', action='write')
  Do ir = 1,nr
     Do jz = 1,nz
        Do kt = 1,nphi-1
           Write(99,'(E18.12)') Bphi(ir,jz,kt)
        End Do
     End Do
  End Do
  Close(99)


  !----------------------------------------------------------
  ! 5. Print timing
  !----------------------------------------------------------
  If (rank == 0) Then
     Write(*,*) "Time spent in total: ", get_elapsed_time(tstart), " seconds"
  End If

  
  Call fin_mpi(.false.)
  
End Program create_bgrid_from_coils
