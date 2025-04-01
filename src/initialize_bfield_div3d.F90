Module initialize_bfield_div3d
  Implicit None
Contains

  Subroutine init_bfield(verbose)
    Use setup_bfield_module
    Use parallel_mod, Only : fin_mpi
    Use bgrid_module, Only : nsym
    Use xdr_routines_mod, Only : nperio
    Implicit None
    Logical, Intent(In) :: verbose

    setup_bfield_verbose = verbose

    Select Case (rmp_type)
    Case ('g')
       Call setup_bfield_g3d
       If (nfp_bfield .ne. 1) Then
          Write(*,*) "Error: bfield method is g but nfp_bfield is not 1",nfp_bfield
          Call fin_mpi(.true.)
       Endif
#ifdef HAVE_FXDR 
    Case ('xdr')
       Call setup_bfield_xdr
       If (nfp_bfield .ne. nperio) Then
          Write(*,*) "Error: nfp_bfield is not equal to nperio for xdr",nfp_bfield,nperio
          Call fin_mpi(.true.)
       Endif
#endif 
    Case ('vmec_coils')
       Call setup_bfield_vmec_coils
    Case ('vmec_coils_to_fil')
       Call setup_bfield_vmec_coils_to_fil
    Case ('bgrid')
       Call setup_bfield_bgrid
       If (nfp_bfield .ne. nsym) Then
          Write(*,*) "Error: nfp_bfield is not equal to nsym for bgrid",nfp_bfield,nsym
          Call fin_mpi(.true.)
       Endif
    Case Default
       If (verbose) Then
          Write(*,*) 'Unknown rmp_type in div3d: ',Trim(rmp_type)
          Write(*,*) 'Current options are:'
          Write(*,*) '''g'''      
          Write(*,*) '''vmec_coils'''
          Write(*,*) '''vmec_coils_to_fil'''
          Write(*,*) '''bgrid'''
#ifdef HAVE_FXDR 
          Write(*,*) '''xdr'''
#endif
       Endif
       Call fin_mpi(.true.)
    End Select



  End Subroutine init_bfield

End Module initialize_bfield_div3d
