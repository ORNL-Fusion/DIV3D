Subroutine find_nearby_triangles_v2(rr,pp,zz,dL,period)

Use kind_mod
Use read_parts_mod
Use math_routines_mod, Only: line_seg_facet_int
Implicit None
Real(real64),Intent(in) :: rr,pp,zz,dL, period
Real(real64) :: xm,ym,zm,R2,P2,Z2,X2,Y2,dd
Integer(int32) :: ipart, itri, icount, itmp
Real(real64), allocatable :: near_part_tmp(:), near_tri_tmp(:)
Integer(int32) :: ii
! test

itmp = sum(ntri_parts(1:nparts))
allocate(near_part_tmp(itmp))
allocate(near_tri_tmp(itmp))

R2 = rr
P2 = pp
Z2 = zz
Do While (pp .lt. 0.d0)
  P2 = P2 + period
Enddo
P2 = mod(P2,period)
X2 = R2*cos(P2)
Y2 = R2*sin(P2)

icount = 1
ic_near = 0
Do ipart=1,nparts  
  Do itri = 1,ntri_parts(ipart)
    xm = xmid(ipart,itri)
    ym = ymid(ipart,itri)
    zm = zmid(ipart,itri)
    
    dd = sqrt( (X2-xm)*(X2-xm) + (Y2-ym)*(Y2-ym) + (Z2-zm)*(Z2-zm) )
    
!    write(*,*) 'this:',dd,dL
!    stop

    if (dd .le. 1.2d0*( dL + dmid(ipart,itri) )) Then    
!    if (dd .le. dL) Then    
      ic_near = ic_near + 1
      near_part_tmp(ic_near) = ipart    
      near_tri_tmp(ic_near) = itri
    endif
    
    icount = icount + 1    
  Enddo ! itri
Enddo ! ipart


allocate(near_part(ic_near),near_tri(ic_near))
near_part = near_part_tmp(1:ic_near)
near_tri = near_tri_tmp(1:ic_near)

deallocate(near_part_tmp,near_tri_tmp)

!write(*,*) 'hjeeree>>>>',ic_near
!stop



End Subroutine find_nearby_triangles_v2
