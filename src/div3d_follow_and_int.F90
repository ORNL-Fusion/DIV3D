!-----------------------------------------------------------------------------
!+ Program to calculate intersections of fieldlines initiated on a flux surface
!  and diffused until they intersect components
!-----------------------------------------------------------------------------
Program div3d_follow_and_int
! Author(s): J.D. Lore - 07/14/2011 - xxx
!
! Modules used:
Use run_settings_namelist, Only : read_run_settings_namelist, ns_line_surf, trace_surface_opt
Use parallel_mod, Only : rank, nprocs, init_mpi, fin_mpi
Use read_parts_mod, Only : read_parts, make_triangles
Use diffusion, Only : diffuse_lines3, diffuse_lines3_worker
Use initialize_bfield_div3d, Only : init_bfield
Use surface_mod, Only : trace_surface
Use initialize_points, Only : init_points_line
Implicit none
Logical :: verbose = .false.
! Variables for timing
Integer :: start = 1, finish = 1, count_rate = 1
Real :: elapsed_time

!----------------------------------------------------------
! 0. Setup
! --Initialize timing
! --Initialize MPI
! --Read namelist
!----------------------------------------------------------
If (rank .eq. 0) Then
   call system_clock(count_rate = count_rate)
   call system_clock(start)
End If

Call init_mpi
If (rank .eq. 0) verbose = .true. 
If (verbose) Write(*,'(/A)') '-------------------------------------------------------------------------'

! Read namelists
Call read_run_settings_namelist(verbose)

!----------------------------------------------------------
! 1. Initialize magnetic field
!----------------------------------------------------------
Call init_bfield(verbose)

!----------------------------------------------------------------
! 2. Load intersection components
!  -- read parts.list file for file names
!  -- read part files
!  -- make trianges from parts
!----------------------------------------------------------------
If (verbose) Write(*,'(A)') ' Reading parts list and part files:'
Call read_parts(verbose)
  
! Make triangles from 2d parts
If (verbose) Write(*,'(/A)') ' Generating 2d part triangles'
Call make_triangles(verbose)

!----------------------------------------------------------------
! Root node traces initial surface and defines starting points
!----------------------------------------------------------------
If (rank .eq. 0) Then

  Write(*,'(A,I0,A)') ' Root node (rank ',rank,') is initializing run'
  Write(*,'(A,I0,A)') ' Total number of processes: ',nprocs

  !----------------------------------------------------------
  ! 3. Trace out a surface (no diffusion)
  !----------------------------------------------------------
  If (trace_surface_opt .eqv. .true.) Then
    Call trace_surface(ns_line_surf)

    !----------------------------------------------------------
    ! 4. Initialize points on surface to carry heat
    !----------------------------------------------------------
    Write(*,*) 'Initializing points along initial surface line'
    Call init_points_line
  Else
    Write(*,*) 'Skipping trace_surface, loading init points from file.'
  Endif
Endif

! Display timing info for prep
If (rank .eq. 0) Then
   call system_clock(finish)
   elapsed_time = Real(finish - start)/Real(count_rate)
   Write(*,*) "Time spent in prep: ", elapsed_time, " seconds"
   call system_clock(start) ! Re-initialize for next call
End If


If (rank .eq. 0) Then

   !----------------------------------------------------------------------------------------
   ! 5. Root node distributes jobs using init points
   !----------------------------------------------------------------------------------------
   
   Write(*,'(/A)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
   write(*,*) 'Beginning MPI line diffusion'
   Call diffuse_lines3
  
Else

   !----------------------------------------------------------------
   ! 5. Additional nodes wait for lines and follow them
   !    Then check for intersections 
   !----------------------------------------------------------------
   Call diffuse_lines3_worker
Endif

! Write timing info
If (rank .eq. 0) Then
   call system_clock(finish)
   elapsed_time = Real(finish - start) / Real(count_rate)
   Write(*,*) "Time spent in following and intersection: ", elapsed_time, " seconds"
End If
 
! Finialize MPI
Call fin_mpi(.false.) ! False means this is a non-error exit

end program div3d_follow_and_int





