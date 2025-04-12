module module_sf_ucm_vf


contains
!===========================================================
subroutine cal_view_factor (                         &
                              hl, wl, rl, ncap,      & !- i
                              svf_r, svf_g, svf_w,   & !- o
                              vf_wwp, vf_wwv,        & !- o
                              vf_gw, vf_wg           & !- o
                            )

implicit none
real(8), intent(in)   :: hl, wl, rl
integer, intent(in)   :: ncap
real(8), intent(out)  :: svf_r, svf_g
real(8), dimension(ncap), intent(out)      :: svf_w
real(8), dimension(ncap), intent(out)      :: vf_gw
real(8), dimension(ncap), intent(out)      :: vf_wg
real(8), dimension(ncap,ncap), intent(out) :: vf_wwp, vf_wwv


svf_r = 1.D0
CALL vf_ground02(                          &
                 hl, wl, rl, ncap,         & !- in
                 svf_g,                    & !- out
                 vf_gw                     & !- out
                 )
CALL vf_wall(                              &
                 hl, wl, rl, ncap,         & !- in
                 svf_w, vf_wg,             & !- out
                 vf_wwp, vf_wwv            & !- out
             )

RETURN
END subroutine cal_view_factor
!==================================================================





!==================================================================
subroutine vf_ground02(                          &
                       hl, wl, rl, ncap,         & !- in
                       svf_g,                    & !- out
                       vf_gw                     & !- out
                       )
implicit none
integer, intent(in)               :: ncap
real(8), intent(in)               :: hl, wl, rl     
real(8), intent(out)              :: svf_g   
real(8), dimension(ncap), intent(out) :: vf_gw   ! ground -> wall
!- local variables
real(8) :: SVF(2), VF(2,ncap)
integer                           :: k, i, n, IS, ID
real(8), dimension(1:ncap+1)      :: z_cap
real(8)                           :: dzcap
real(8)                           :: x1, x2, Y1, Y2
real(8)  :: PI
integer  :: ILAT, ILON, NLAT, NLON
real(8)  :: DLAT, DLON, LAT, LON

integer, parameter :: dim_num = 3
integer, parameter :: REF_SURFACE = 8
real(8) :: F, G, H
real(8), dimension(1:REF_SURFACE*2)          :: AA, BB, CC, DD
real(8), dimension(REF_SURFACE*2,DIM_NUM)  :: PS
logical, dimension(REF_SURFACE*2)          :: intersect
real(8), dimension(DIM_NUM)                  :: P, P1, P2, PCOR
real(8)    :: x0, Y0, Z0, phanbiet, DIRS(2,4), UU
integer    ::    inDS, IX, IY, NRW, IGP
real(8), dimension(:,:), allocatable   :: SIGNAL, SIGNAL_DIR
LOGICAL, dimension(:,:), allocatable   :: MASK
real(8), dimension(:), allocatable     :: XXS, XXE, YYS, YYE
real(8) :: AREA1, AREA2
CHARACTER(LEN=10), dimension(2)        :: FILENAME
!--------------------------------------------------
!    ________       ________
!    |       |      |       |
!    |       |      |       |
!    |       |      |       |
!    |_______|      |_______|
!
!        1      2
!    ________       ________
!    |       |      |       |
!    |       |  1   |       |
!    |       |      |       |
!    |_______|      |_______|
!--------------------------------
PI =  4.d0*datan(1.d0)
n = ncap
DZCAP = hl/dfloat(ncap)
z_cap(1) = 0.0d0
do k = 2, ncap+1
   z_cap(k) = z_cap(k-1) + dzcap
end do


FILENAME(:) = (/"VITRI1","VITRI2"/)
do IGP= 1, 2
   OPEN(100+IGP,FILE=".GROUND_VIEW_FACTOR_"//TRIM(FILENAME(IGP)) )

   UU = RL+WL
   NRW = 10
   allocate( XXS( NRW ) )
   allocate( XXE( NRW ) )
   allocate( YYS( NRW ) )
   allocate( YYE( NRW ) )

   XXS(1)  = -(DFLOAT(NRW/2 ) * UU  -  RL/2.D0)
   XXE(1)  = XXS(1) + WL
   do I = 2, NRW
      XXS(I) = XXS(I-1) + UU
      XXE(I) = XXS(I) + WL
   enddo

   YYS(:)  = XXS(:)
   YYE(:)  = XXE(:)

   do K = 1, REF_SURFACE*2
      AA (K) = DFLOAT( MOD((K+1)/2,2) )       
      BB (K) = DFLOAT( MOD((K-1)/2,2) )
      CC (K) = 0.D0                           
      DD (K) = -(RL/2.D0 + UU * DFLOAT((K-1)/4))  &
                * DCOS(DFLOAT(K) * PI) 
   enddo

   IF(IGP==2) then
       IF(.TRUE.) then

         YYS(1)  = -(DFLOAT(NRW/2 ) * UU  +  WL/2.D0 )
         YYE(1)  = YYS(1) + WL
         do I = 2, NRW
            YYS(I) = YYS(I-1) + UU
            YYE(I) = YYS(I) + WL

         enddo

         !-----Ground 1
         do I = 1, REF_SURFACE*2
            IF(MOD(I,4)==0) then
               do K = I-1, I 
                  DD (K) = -(WL/2.D0 + RL + UU  &
                           * DFLOAT((K-1)/4))   &
                           * DCOS(DFLOAT(K) * PI) 

                enddo
            END IF
         enddo
      END IF
   END IF


   IF(.FALSE.) then
      XXS(1)  = -( DFLOAT(NRW/2 ) * UU -RL / 2.D0 + 4.D0)
      XXE(1)  = XXS(1) + WL
      do I = 2, NRW
         XXS(I) = XXS(I-1) + UU
         XXE(I) = XXS(I) + WL
      enddo

      !-----Ground 1
      do I = 1, REF_SURFACE*2
         IF(MOD(I-1,4)==0) then
            do K = I, I+1 
               DD (K) = DD(K) - 4.D0
            enddo
         END IF
      enddo
   END IF



   DLAT = 3.D0
   DLON = 3.D0

   NLAT =  inT(90.D0/DLAT) + 1
   NLON =  inT(360.D0/DLON) 
   allocate( SIGNAL     (1:NLON,1:NLAT) )
   allocate( SIGNAL_DIR (1:NLON,1:NLAT) )
   allocate( MASK       (1:NLON,1:NLAT) )
   SIGNAL(:,:)      = -1.D0
   SIGNAL_DIR(:,:)  = -1.D0


   do ILAT = 1, NLAT
      LAT       =  DFLOAT(ILAT-1)*DLAT*PI/180.D0
      do ILON = 1, NLON
         LON       =  DFLOAT(ILON-1)*DLON*PI/180.D0 + DLON*PI/180.D0/2.D0
         Z0 = DSin(LAT) 
         X0 = DCOS(LON) * DCOS(LAT)
         Y0 = DSin(LON) * DCOS(LAT)

         P1(1:3) = (/0.D0, 0.D0, 0.D0 /)
         P2(1:3) = (/X0, Y0, Z0 /)
          
         PCOR(1:3) = (/0.D0, 0.D0,0.D0/)

         P1(1:3)   = P1(1:3) + PCOR(1:3)
         P2(1:3)   = P2(1:3) + PCOR(1:3)


         do IS = 1, REF_SURFACE*2, 1   ! SURFACE

            IF(SIGNAL(ILON,ILAT) .NE. -1.D0) EXIT

            CALL line_exp2par_3d ( p1, p2, f, g, h, x0, y0, z0 )

            CALL plane_imp_line_par_int_3d (                    &
                         AA(IS), BB(IS), CC(IS), DD(IS),        &
                                             x0, y0, z0,        &
                                             f, g, h,           &
                                             intersect(IS),     &
                                             PS(IS,:) )

            IF(inTERSECT(IS) .and. doT_PRODUCT(P2,PS(IS,:)) > 0.D0) then
               do K = 2, ncap+1
                  IF( PS(IS,3) < Z_CAP(K) .and. PS(IS,3) >= Z_CAP(K-1)) then

                     ! vertical wall
                     IF(AA(IS) == 1.D0) then
                        do I = 1,NRW
                           IF(PS(IS,2) >= YYS(I) .and. PS(IS,2) <= YYE(I) ) then
                              SIGNAL(ILON,ILAT) = DFLOAT(K-1)
                           END IF
                        enddo
                     END IF  
            
                     ! horizontal walls
                     IF(BB(IS) == 1.D0) then
                        do I = 1,NRW
                           IF(PS(IS,1) >= XXS(I) .and. PS(IS,1) <= XXE(I) ) then
                              SIGNAL(ILON,ILAT) = DFLOAT(K-1)
                           END IF
                        enddo
                     END IF
                  END IF
               enddo  ! enddo K
          END IF
      enddo ! enddo DIRECTION ID
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !write (500+IGP, '(2x,4F10.3,2x)' ) LAT*180.D0/PI, LON*180.D0/PI,  &
      !                             SIGNAL(ILON,ILAT), SIGNAL_DIR(ILON,ILAT)
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   enddo
enddo


write(20,*) '--- ground view factor ------------------------- '
MASK = SIGNAL(:,:) .NE. -1.D0 
SVF(IGP)  = 1.d0 - DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))
write(20,'(A, F7.3)')  '    SKY VIEW:    ', 1.d0 - DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))

do K = 1, ncap
   MASK = SIGNAL(:,:) .EQ.  DFLOAT(K) 
   write(20,'(A, F7.3)') '    VF -> WALL  :  ', DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))
   VF(IGP,K) = DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))
enddo

deallocate(MASK, SIGNAL, SIGNAL_DIR)
deallocate(XXS, XXE, YYS, YYE)

enddo

   write(20,*) '--- ground view factor average------------------- '
     AREA1 = RL*RL
     AREA2 = RL*WL
     SVF_G  = (SVF(1)*AREA1 + 2.D0*SVF(2)*AREA2) / (AREA1+2.D0*AREA2)
   write(20,'(A, F7.3)')  '    SKY VIEW:    ', SVF_G


   do K = 1, ncap
    VF_GW(K)  = (VF(1,K)*AREA1 + 2.D0*VF(2,K)*AREA2) / (AREA1+2.D0*AREA2)
    write(20,'(A, F7.3)') '    VF -> WALL  :  ', VF_GW(K)
   enddo

  
   write(20, '(A, F7.3)') '  CHECK:  ', SVF_G + SUM(VF_GW)

  return
  end subroutine vf_ground02
!===================================================================











!===================================================
subroutine vf_wall(                              &
                       hl, wl, rl, ncap,         & !- in
                       svf_w, vf_wg,             & !- out
                       vf_wwp, vf_wwv            & !- out
                   )

implicit none
integer, intent(in)                 :: ncap
real(8), intent(in)                 :: hl, wl, rl    
real(8), dimension(1:ncap),intent(out)   :: svf_w  ! wall (layers)
real(8), dimension(1:ncap),intent(out)   :: vf_wg  ! wall -> ground
real(8), dimension(1:ncap,1:ncap),intent(out) :: vf_wwp, vf_wwv 

!- local variables
real(8), dimension(1:ncap+1)        :: z_cap
real(8)                             :: dzcap
real(8)  :: PI
real(8)  :: DLAT, DLON, LAT, LON
integer, parameter :: dim_num = 3
integer, parameter :: REF_SURFACE = 8
  
real(8) :: A, B, C, D, F, G, H
real(8), dimension(1:REF_SURFACE)          :: AA, BB, CC, DD
real(8), dimension(REF_SURFACE,DIM_NUM)    :: PS
logical, dimension(REF_SURFACE)            :: intersect
real(8), dimension(DIM_NUM)                  :: p, P1, P2
real(8)    :: x0, Y0, Z0   
integer    ::    ik, IW, ILAT, ILON, k, i, n, IS, NLAT, NLON
real(8), dimension(:,:), allocatable   :: SIGNAL, SIGNAL_DIR
logical, dimension(:,:), allocatable   :: MASK


PI =  4.d0*datan(1.d0)
N = ncap
dzcap = HL / dfloat(ncap)
z_cap(1) = 0.0d0

do k = 2, ncap + 1
   z_cap(k) = z_cap(k-1) + dzcap
end do

DLAT = 3.D0
DLON = 3.D0

do IK = 1, N
   P1(1:3) = (/0.D0, 0.D0, 0.D0 /)
   NLAT =  inT(90.D0/DLAT) + 1
   NLON =  inT(360.D0/DLON) 
   allocate(SIGNAL(1:NLON,1:NLAT))
   allocate(SIGNAL_DIR(1:NLON,1:NLAT))
   allocate(MASK(1:NLON,1:NLAT))

   SIGNAL(:,:)      = -1.D0
   SIGNAL_DIR(:,:)  = -1.D0

   do ILAT = 1, NLAT
      LAT       =  DFLOAT(ILAT-1)*DLAT*PI/180.D0
      do ILON = 1, NLON
         LON       =  DFLOAT(ILON-1)*DLON*PI/180.D0
         Z0 = DSin(LAT) 
         X0 = DCOS(LON) * DCOS(LAT)
         Y0 = DSin(LON) * DCOS(LAT)
         P2(1:3) = (/X0, Y0, Z0 /)

         ! WALL SURFACE
         AA(1)   =  0.D0
         BB(1)   =  0.D0
         CC(1)   =  1.D0
         DD(1)   =  -RL

         ! MAT PHANG SONG SONG VOI NORMAL VECTOR CUA WALL POinT
         AA(2:5)   = 1.D0
         BB(2:5)   = 0.D0
         CC(2:5)   = 0.D0
         DD(2:5)   = (/ -WL -RL - WL /2.D0 - RL, - WL /2.D0  - RL, &
                           WL /2.D0 +RL, WL + RL + WL/2.D0 +RL/)  

         ! MAT PHANG VUONG GOC VOI NORMAL VECTOR CUA WALL POinT
         AA(6)   = 0.D0
         BB(6)   = 0.D0
         CC(6)   = 1.D0
         DD(6)   = -(RL+WL+RL)

         ! MAT PHANG ROAD
         AA(7)    = 0.D0
         BB(7)    = 1.D0
         CC(7)    = 0.D0
         DD(7)    = -Z_CAP(IK) - DZCAP/2.D0

         do IS = 7, 1,-1 ! SURFACE
            a = AA(IS)
            b = BB(IS)
            c = CC(IS)
            d = DD(IS) 
         
            CALL line_exp2par_3d ( p1, p2, f, g, h, x0, y0, z0 )
            CALL plane_imp_line_par_int_3d ( a, b, c, d,        &
                                             x0, y0, z0,        &
                                             f, g, h,           &
                                             intersect(IS),     &
                                             PS(IS,:) )

            IF(inTERSECT(IS)) then  !
            IF(doT_PRODUCT(P2,PS(IS,:)) > 0.D0) then

            !----view to road---------------------
            IF(IS==7) then   ! ground-> signal 1
                       SIGNAL(ILON,ILAT)     = 0.D0
                       SIGNAL_DIR(ILON,ILAT) = 0.D0
               END IF


               IF(IS == 6) then  ! for signal 2: perpendicular surface
                  do K = 1, ncap   ! canopy layer number
 
                     IF( PS(IS,2)<=Z_CAP(IK)+DZCAP/2.D0-Z_CAP(K).and.  &
                         PS(IS,2)>=Z_CAP(IK)+DZCAP/2.D0-Z_CAP(K)-DZCAP) then 
                       
                         do IW = -3, 3
                            if( PS(IS,1)<=DFLOAT(IW)*(RL+WL)+WL/2.D0 .and.&
                                PS(IS,1)>=DFLOAT(IW)*(RL+WL)-WL/2.D0) then
                               SIGNAL(ILON,ILAT) = DFLOAT(K) 
                               SIGNAL_DIR(ILON,ILAT) = 1.D0
                           endif
                        enddo
                     endif
                  enddo
               END IF


               IF(IS>=2 .and. IS <= 5) then  ! FOR signal 3: parallel surface
                  do K = 1, ncap 
                     IF(PS(IS,2)<=Z_CAP(IK)+DZCAP/2.D0-Z_CAP(K) .and.   &
                        PS(IS,2)>=Z_CAP(IK)+DZCAP/2.D0-Z_CAP(K)-DZCAP) then
                     
                        IF(PS(IS,3) >= RL .and. PS(IS,3) <= (RL+WL ) ) then
                            SIGNAL(ILON,ILAT)     = DFLOAT(K)
                            SIGNAL_DIR(ILON,ILAT) = 2.D0
                        endif
                     endif
                  enddo
               endif

               if(IS==1) then  ! surface perpendicular to normal vector of wall point

                 do K = 1, ncap ! for signal 2
                   if( PS(IS,2)<=Z_CAP(IK) + DZCAP/2.D0-Z_CAP( K ) .and.        &
                       PS(IS,2)>=Z_CAP(IK) + DZCAP/2.D0-Z_CAP( K ) - DZCAP ) then 
                         do IW = -5, 5
                            IF( PS(IS,1)  <=  DFLOAT(IW)*(RL+WL) + WL/2.D0 .and.     &
                                PS(IS,1)  >=  DFLOAT(IW)*(RL+WL) - WL/2.D0 ) then
                                SIGNAL(ILON,ILAT) = DFLOAT(K) 
                                SIGNAL_DIR(ILON,ILAT) = 1.D0

                            END IF
                         enddo
                     endif
                  enddo
               endif
            endif
         endif
      enddo
   enddo
enddo

write(20,*) 'LAYER------------------------- ',   IK
MASK = SIGNAL(:,:) .NE. -1.D0
SVF_W( IK )  =  1.d0 - DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))
write(20,'(A, F7.3)')  '    SKY VIEW:    ', SVF_W(IK)

MASK = SIGNAL(:,:) .EQ.  0.D0 
VF_WG( IK )  = DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))
MASK = SIGNAL(:,:) .EQ.  0.D0 
write(20,'(A, F7.3)')  '    GROUND VIEW: ', VF_WG( IK )
       
do K = 1, ncap

   MASK = SIGNAL(:,:) .EQ.  DFLOAT(K) .and. SIGNAL_DIR(:,:) .EQ. 1.D0
   VF_WWP( IK,K )  = DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))

   MASK = SIGNAL(:,:) .EQ.  DFLOAT(K) .and. SIGNAL_DIR(:,:) .EQ. 2.D0
   VF_WWV( IK,K )  =  DFLOAT(COUNT(MASK)) / DFLOAT(SIZE(MASK))
   write(20,'(A, 2F10.4)') '    VF -> P -> V       :  ',VF_WWP( IK,K ), VF_WWV( IK,K )
enddo
   
write(20,'(A, F7.3)') '    CHECK       :  ',  SUM( VF_WWV( IK,: ) ) &
                                            + SUM( VF_WWP( IK,: ) ) &
                                            + VF_WG( IK ) + SVF_W(IK)

deallocate(MASK, SIGNAL, SIGNAL_DIR)   
enddo


return
end subroutine vf_wall
!===================================================================











subroutine plane_imp_line_par_int_3d ( a, b, c, d, x0, y0, z0, f, g, h, &
                                       intersect, p )
!************************************************************************
!
!! PLANE_IMP_LinE_PAR_inT_3D: intersection ( impl plane, param line ) 
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F*F + G*G + H*H = 1,
!    and F nonnegative.
!
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmers Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420,
!    page 111.
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the implicit plane parameters.
!
!    Input, real ( kind = 8 ) X0, Y0, Z0, F, G, H, parameters that define the
!    parametric line.
!
!    Output, logical ( kind = 4 ) inTERSECT, is TRUE if the line and the plane
!    intersect.
!
!    Output, real ( kind = 8 ) P(3), is a point of intersection of the line
!    and the plane, if inTERSECT is TRUE.
!
implicit none
integer, parameter    :: dim_num = 3
real(8), intent(in)   :: a, b, c, d, x0, y0, z0, f, g, h
real(8), intent(out)  :: p(dim_num)
logical, intent(out) :: intersect

real ( 8 ) denom
real ( 8 ) norm1
real ( 8 ) norm2
real ( 8 ) t
real ( 8 ), parameter :: tol = 0.00001D+00
 
!
!  Check.
!
  norm1 = sqrt ( a * a + b * b + c * c )

  if ( norm1 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LinE_PAR_inT_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop 1
  end if

  norm2 = sqrt ( f * f + g * g + h * h )

  if ( norm2 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LinE_PAR_inT_3D - Fatal error!'
    write ( *, '(a)' ) '  The line direction vector is null.'
    stop 1
  end if

  denom = a * f + b * g + c * h
!
!  The line and the plane may be parallel.
!
  if ( abs ( denom ) < tol * norm1 * norm2 ) then

    if ( a * x0 + b * y0 + c * z0 + d == 0.0D+00 ) then
      intersect = .true.
      p(1) = x0
      p(2) = y0
      p(3) = z0
    else
      intersect = .false.
      p(1:dim_num) = 0.0D+00
    end if
!
!  If they are not parallel, they must intersect.
!
  else

    intersect = .true.
    t = - ( a * x0 + b * y0 + c * z0 + d ) / denom
    p(1) = x0 + t * f
    p(2) = y0 + t * g
    p(3) = z0 + t * h

  end if

  return
end subroutine





subroutine line_exp2par_3d ( p1, p2, f, g, h, x0, y0, z0 )

!*****************************************************************************80
!
!! LinE_EXP2PAR_3D converts a line from explicit to parametric form in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F^2 + G^2 + H^2 = 1, and F nonnegative.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(3), P2(3), two points on the line.
!
!    Output, real ( kind = 8 ) F, G, H, X0, Y0, Z0, the parametric parameters
!    of the line.
!
  implicit none

  integer, parameter   :: dim_num = 3
  real(8), intent(in)  :: p1(dim_num), p2(dim_num)
  real(8), intent(out) :: x0, y0, z0
  real(8), intent(out) :: F, G, H

  logical              :: line_exp_is_degenerate_nd
  real(8)              :: norm


  
  line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

   if ( line_exp_is_degenerate_nd ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'LinE_EXP2PAR_3D - Warning!'
     write ( *, '(a)' ) '  The line is degenerate.'
   end if

  x0 = p1(1)
  y0 = p1(2)
  z0 = p1(3)

  f = p2(1) - p1(1)
  g = p2(2) - p1(2)
  h = p2(3) - p1(3)

  norm = sqrt ( f * f + g * g + h * h )

  if ( norm /= 0.D0 ) then
    f = f / norm
    g = g / norm
    h = h / norm
  end if

  if ( f < 0.D0 ) then
    f = -f
    g = -g
    h = -h
  end if

  return
end subroutine line_exp2par_3d




END MODULE module_sf_ucm_vf
