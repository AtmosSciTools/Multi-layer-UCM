module module_initialize

use module_params
use module_urb_ini
use module_vars, only:                           &
         ZR, ROOF_WIDTH, ROAD_WITH, AH, FRC_URB, CAPR, CAPB, CAPG,     &
         AKSR, AKSB, AKSG, ALBR, ALBB, ALBG, Z0B, Z0G, Z0R, TRLEND,    &
         TBLEND, TGLEND, DDZR, DDZB, DDZG,                             &
         THERMAL_INSOL_ROOF, THERMAL_INSOL_WALL, RLNU, BLNU, GLNU,     &
         BOUNDR, BOUNDB, BOUNDG, AHOPTION, FRCURBOPT,                  &
         AHDIUPRF,                                                     &
         MOLS, MOLR, MOLB, MOLG, MOLC,     &   ! Monin-Obukhov length
         obs_vars

use module_sf_ucm_vf

contains

!==============================================================================
 subroutine initialize(                                                &
                          s0, elevation,                               & !- out
                          adi, alpha_s, alpha_c,                       & !- out
                          km, kh,                                      & !- out
                          tau_u, tau_v, tau_t,                         & !- out
                          tau_qv,                                      & !- out
                          z, z_t, dz, dz_t,                            & !- out
                          svf_r, svf_g, svf_w,                         & ! -out
                          vf_gw, vf_wg, vf_wwp, vf_wwv,                & !- out
                          sw, sd, ss, lw,                              & !- out
                          sw_ra_heating, lw_ra_heating,                & !- out
                          kce, canop_a, canop_m, q, qa, qev,           & !- out
                          year, month, day, hour, minute, second,      & !- out
                          tssl, tgsl, ts_g, ts_r, tw_r, tw_g, ts, tw,  & !- out
                          hl, wl, rl,                                  & !- out
                          u, v, p, t,                                  & !- out
                          uu, vv, tt,                                  & !- out
                          qv, rho_a, rh                                & !- out
                      )
!==============================================================================

implicit none

real(8), intent(out)                         :: s0, elevation
real(8), intent(out)                         :: adi, alpha_s, alpha_c
real(8), dimension(kms:kme), intent(out)     :: km, kh
real(8), intent(out)                         :: tau_u, tau_v, tau_t
real(8), intent(out)                         :: tau_qv
real(8), intent(out)                         :: sw, sd, ss, lw
real(8), dimension(kms:kme), intent(out)     :: sw_ra_heating, lw_ra_heating
integer, intent(out)                         :: year, month, day,      &
                                                  hour, minute, second
real(8), dimension(kms:kme+1), intent(out)   :: z
real(8), dimension(kms:kme), intent(out)     :: z_t, dz, dz_t
real(8), dimension(kms:kme), intent(out)     :: u, v, t, qv
real(8), dimension(kms:kme), intent(out)     :: p

real(8), dimension(kms:kme), intent(out)     :: uu, vv, tt
real(8), dimension(kms:kme), intent(out)     :: rho_a, rh
! urban canopy group
integer, intent(out)                            :: kce
real(8), intent(out)                            :: hl, wl, rl
real(8), dimension(kms:kme), intent(out)        :: canop_a, canop_m
real(8), intent(out)                            :: svf_r, svf_g
real(8), intent(out)                            :: ts_r, ts_g, tssl
real(8), dimension(1:nw), intent(out)           :: tgsl
real(8), dimension(kms:kme), intent(out)        :: q, qa, qev
real(8), dimension(:), allocatable, intent(out)      :: tw_r, tw_g
real(8), dimension(:), allocatable, intent(out)      :: svf_w
real(8), dimension(:), allocatable, intent(out)      :: vf_gw
real(8), dimension(:), allocatable, intent(out)      :: vf_wg
real(8), dimension(:,:), allocatable, intent(out)    :: vf_wwp, vf_wwv
real(8), dimension(:,:), allocatable, intent(out)    :: ts
real(8), dimension(:,:,:), allocatable, intent(out)  :: tw


call vars_zero_clear (                                  &
                      s0, elevation,                    & !- out
                      adi, alpha_s, alpha_c,            & !- out
                      z, z_t, dz, dz_t,                 & !- out
                      u, v, p, t,                       & !- out
                      q, qa, qev,                       & !- out
                      uu, vv, tt,                       & !- out
                      qv, rho_a,                        & !- out
                      km, kh,                           & !- out
                      tau_u, tau_v, tau_t,              & !- out
                      tau_qv,                           & !- out
                      sw, sd, ss, lw,                   & !- out
                      sw_ra_heating, lw_ra_heating      & !- out
                      )
year   = start_year
month  = start_month
day    = start_day
hour   = start_hour
minute = start_minute
second = start_second

call read_urbparm ()

call assign_urb_parm ( UTYPE )

call grid( z, z_t, dz, dz_t )              !- set grid

call dyn_init( z_t, u, v, t, rho_a, p )

call moist_init(                                &
                       rho_a, z, p, t,          & !- in
                       rh, qv, qev              & !- out
                )

call read_landuse (sltype)
!..........................
!-for urban condition

call build( z, hl, wl, rl, kce )

call ucm_vars_allocate (                         &
                       kce, svf_w, vf_gw,        &
                       vf_wg,vf_wwp, vf_wwv,     &
                       ts, tw, tw_r, tw_g        &
                        )

call canop (                            &
              z, dz, hl, wl, rl, kce,   & !- in
              canop_a, canop_m          & !- inout
            )

!- calculate the types of view factor, included
!  sky view factor of each surface at each layer
!  view factor of each surface at each layer to other surface
!  by use building, and road feature
!  independent to time step

call cal_view_factor (                            &
                        hl, wl, rl, kce-2,        & !- in
                        svf_r, svf_g, svf_w,      & !- out
                        vf_wwp, vf_wwv,           & !- out
                        vf_gw, vf_wg              & !- out
                      )

vf_gw(:) = vf_gw(:) / 4.D0
vf_wwv(:,:) = vf_wwv(:,:) / 2.D0

call ts_all ( kce, tssl, tgsl, ts_r, ts_g, tw_r, tw_g, ts, tw  )  ! - out

return
end subroutine initialize
!=========================================================================







!======================================================================
subroutine vars_zero_clear (                                  &
                              s0, elevation,                  & !- out
                              adi, alpha_s, alpha_c,          & !- out
                              z, z_t, dz, dz_t,               & !- out
                              u, v, p, t,                     & !- out
                              q, qa, qev,                     & !- out
                              uu, vv, tt,                     & !- out
                              qv, rho_a,                      & !- out
                              km, kh,                         & !- out
                              tau_u, tau_v, tau_t,            & !- out
                              tau_qv,                         & !- out
                              sw, sd, ss, lw,                 & !- out
                              sw_ra_heating, lw_ra_heating    & !- out
                             )
implicit none
real(8), intent(out)       :: s0, elevation
real(8), intent(out)       :: adi, alpha_s, alpha_c
real(8), intent(out)       :: sw, sd, ss, lw
real(8), intent(out)       ::  tau_u, tau_v, tau_t, tau_qv
real(8), dimension(kms:kme+1), intent(out) :: z
real(8), dimension(kms:kme), intent(out)   :: z_t, dz, dz_t
real(8), dimension(kms:kme), intent(out)   :: u, v
real(8), dimension(kms:kme), intent(out)   :: p
real(8), dimension(kms:kme), intent(out)   :: t, qv
real(8), dimension(kms:kme), intent(out)   :: uu, vv, tt
real(8), dimension(kms:kme), intent(out)   :: rho_a
real(8), dimension(kms:kme), intent(out)   :: km, kh
real(8), dimension(kms:kme), intent(out):: sw_ra_heating,lw_ra_heating
real(8), dimension(kms:kme), intent(out)   :: q, qa, qev

q(:)   = 0.d0
qa(:)  = 0.d0
qev(:) = 0.d0

z(:)   = 0.d0
z_t(:) = 0.d0
dz(:)  = 0.d0
dz_t(:) = 0.d0

u(:) = 0.d0
v(:) = 0.d0
p(:) = 0.d0
t(:) = 0.d0
qv(:) = 0.d0
rho_a(:) = 0.d0

uu(:) = 0.0d0
vv(:) = 0.0d0
tt (:) = 0.0d0


km(:) = 0.d0
kh(:) = 0.d0

tau_u  = 0.d0
tau_v  = 0.d0
tau_t  = 0.d0
tau_qv = 0.d0

!- variable for radiation model

sw   = 0.d0
sd   = 0.d0
ss   = 0.d0
lw   = 0.d0
sw_ra_heating(:) = 0.d0
lw_ra_heating(:) = 0.d0

s0        = 0.d0
elevation = 0.d0
adi       = 0.d0
alpha_s   = 0.d0
alpha_c   = 0.d0
MOLS      = 0.d0
MOLR      = 0.d0
MOLB      = 0.d0
MOLG      = 0.d0
MOLC      = 0.d0

end subroutine vars_zero_clear
!==============================================================







!===============================================
subroutine grid( z, z_t, dz, dz_t )  !- inout
implicit none
real(8),dimension(kms:kme+1),intent(inout) :: z
real(8),dimension(kms:kme),intent(inout) :: z_t, dz, dz_t

!- local variable
integer    :: k
real(8)    :: dza
integer    :: find, iostatus, nrow
!++++++++++++++++++++++++++++++++++++++++++++++++

select case( ow_nesting )
case(1)

   !++++++++++++++++++++++++
   print*, 'Read grid from data'
   print*, trim(ics_path)
   find = 99
   open ( find, file=trim(ics_path)//"/grid.dat",          &
          access='SEQUENTIAL',                             &
          status='OLD', action='READ', position='REWIND',  &
          iostat=IOSTATUS )

   if(IOSTATUS>0) stop 'ERROR  open grid file'

   nrow = 0
   do
      read(find,'()',end=100)
      nrow = nrow + 1
   end do
   100 close(find)
   if (nrow /= kme) stop "Number of grid in namelist is wrong"

   open(find, file=trim(ics_path)//"/grid.dat")

   do k = 1,kme
      read(find,*) z_t(k)
   end do
   close(find)
   !++++++++++++++++++++++++


case(0)

   select case( grid_type )
   case(1)    ! case of unequal grid

      do k = kms,kme
         dz_t(k) = exp(float(k)/float(kme)*4.d0)
      end do
      dza = z_top / sum(dz_t(kms:kme-1))
      dz_t(:) = dz_t(:) * dza


   case(0)        ! case of equal grids

      do k = kms,kme
         dz_t(k) = z_top / float(kme-1)
      end do


   case default

      write(6,*)"erorr ! check 'grid_type' "
      stop

   end select

   z_t(1) = -dz_t(1)  / 2.d0
   z_t(2) =  dz_t(1)  / 2.d0

   do k = 3, kme
      z_t(k) = z_t(k-1) + dz_t(k-1)
   end do

end select

z(2:kme) = (z_t(1:kme-1) + z_t(2:kme))/2.
z(1) = z_t(1)*2
z(kme+1) = z_t(kme) + z_t(kme) - z_t(kme-1)

dz(kms:kme)     = z(kms+1:kme+1) - z(kms:kme)
dz_t(kms:kme-1) = z_t(kms+1:kme) - z_t(kms:kme-1)
dz_t(kme) = dz_t(kme-1)

!do k = kms, kme
!   print*, z_t(k), z(k), dz(k), dz_t(k)
!end do
return
end subroutine grid
!======================================================================










!=================================================================
subroutine dyn_init( z_t, u, v, t, rho_a, p )

implicit none
!- inout variables
real(8), dimension(kms:kme), intent(in)      :: z_t
real(8), dimension(kms:kme), intent(out)   :: u, v, t,rho_a, p
!- local variables
integer ::  k

do k = kms, kme
   u(k) = Ug
   v(k) = Vg
   t(k) = 0.d0
   rho_a(k) = 1.2d0
   p(k) = pb + ( - rho_a(k) * g * z_t(k) ) * 10.d0**(-2)   ![hPa]
end do

return
end subroutine dyn_init
!=================================================================













!====================================================================
subroutine moist_init(                      &
                       rho_a, z, p, t,      & !- in
                       rh, qv, qev          & !- out
                      )
implicit none

real(8), dimension(kms:kme+1), intent(in)            :: z
real(8), dimension(kms:kme), intent(in)              :: rho_a
real(8), dimension(kms:kme), intent(in)              :: p, t
real(8), dimension(kms:kme), intent(out)             :: rh
real(8), dimension(kms:kme), intent(inout)           :: qv, qev
!- local variables
integer    :: k
real(8)    :: pres    !- pressure [hPa]
real(8)    :: e_sat   !- saturated water vapor pressure [hPa]
real(8)    :: qv_sat  !- saturated mixing ratio of water vapor[kg/kg]
real(8)    :: temp, theta
real(8)    :: xx1, xx2
real(8), dimension(kms:kme)   :: zc
!-- local parameters
real(8), parameter      :: aa = 7.5d0
real(8), parameter      :: bb = 237.3d0    !-- for water sureface
real(8) :: sumqv


qev(kms:kme)  = 0.0d0  ! initialize evaporation source by each layer

do k = kms, kme

   zc(k) = 0.5d0 * ( z(k) + z(k+1) )
   !--------------------------
   ! Set rh
   if ( zc(k) < 1500.d0 ) then
      if ( zc(k) < 0.d0 ) then
         rh(k) = 0.d0
      else
         rh(k) = 80.d0 + ( 60.d0 -  80.d0 ) * zc(k) / ( 1500.d0 - 0.d0 )
        !rh(k) = 50.d0 + ( 20.d0 -  50.d0 ) * zc(k) / ( 1500.d0 - 0.d0 )
      end if
   else
     !rh(k) = 60.d0
      rh(k) = 20.d0

   endif

   !--------------------------
   ! convert to qv from rh
   theta  = tb + ganma * zc(k) + t(k)                      ![K]
   temp   = theta - 0.01d0 * zc(k)
   pres   = p(k)                                           ![hPa]
   !e_sat  = 6.1078d0*10.d0**(aa*(temp-273.15d0)/( bb +(temp - 273.15d0)))
   e_sat  = 6.11d0*10.d0**(7.5d0*(temp-273.15)/((temp-273.15)+237.3d0))
   qv_sat = 0.622d0 * ( e_sat / pres) / ( 1.d0 - 0.378d0 * (e_sat / pres))
   !* 1000 [kg/kg]
   qv(k)  = rh(k)/100.d0 * qv_sat                          !* 4.d0 ![kg/kg]

enddo

!-- for check
!sumqv = 0.d0
!do k = kms, kme
!   sumqv = sumqv + qv(k)
!end do
return
end subroutine moist_init
!==============================================================================






!==================================================
subroutine build( z, hl, wl, rl, kec )  !- inout

implicit none

real(8), dimension(kms:kme+1), intent(in)  :: z
real(8), intent(inout)   :: hl, wl, rl
integer, intent(inout)   :: kec
integer                  :: k

hl =   ZR_TBL         (UTYPE)
wl =   ROOF_WIDTH_TBL (UTYPE)
rl =   ROAD_WITH_TBL  (UTYPE)
!print*, hl, wl, rl
!print*, frc_urb, (wl+rl)/frc_urb
! calculate actual canyon width
rl = (wl+rl)/frc_urb - wl
!print*, hl, wl, rl

kec = 3
do k=kme, kms,-1
   if( z(k) <= hl ) exit
   kec = k
end do

write(20,'(A,F7.2)') 'BUILDING HEIGHT----', HL
write(20,'(A,F7.2)') 'BUILDING WIDTH-----', WL
write(20,'(A,F7.2)') 'ROAD WIDTH     ----', RL
print*, "number of uc layer", kec-2
write(20,*) '   number of urban canopy layers:', kec-2
write(20,'()')

if(kec-2<1) stop 'ERROR: number of UC layers must be > 1'
return
end subroutine build
!====================================================================








!=========================================================
subroutine ucm_vars_allocate ( kce,                 &  !-i
                svf_w, vf_gw, vf_wg,vf_wwp, vf_wwv, &  !-o
                ts, tw, tw_r, tw_g  )

implicit none
integer, intent(in)                                 :: kce
real(8), dimension(:), allocatable, intent(out)     :: tw_r, tw_g
real(8), dimension(:), allocatable, intent(out)     :: svf_w
real(8), dimension(:), allocatable, intent(out)     :: vf_gw
real(8), dimension(:), allocatable, intent(out)     :: vf_wg
real(8), dimension(:,:), allocatable, intent(out):: vf_wwp, vf_wwv
real(8), dimension(:,:), allocatable, intent(out)   :: ts
real(8), dimension(:,:,:), allocatable, intent(out) :: tw

allocate( tw_r( 1:RLNU ) )
allocate( tw_g( 1:GLNU ) )
allocate( ts(1:4, 1:kce-2) )
allocate( tw(1:4, 1:kce-2, 1:BLNU) )
allocate( svf_w(1:kce-2) )
allocate( vf_gw(1:kce-2) )
allocate( vf_wg(1:kce-2) )
allocate( vf_wwp(1:kce-2,1:kce-2) )
allocate( vf_wwv(1:kce-2,1:kce-2) )

return
end subroutine
!============================================================






!============================================================
subroutine canop(                          &
                  z, dz, hl, wl, rl, kce,  & !- in
                  canop_a, canop_m         & !- inout
                )
implicit none
integer, intent(in)                        :: kce
real(8), dimension(kms:kme), intent(in)    :: dz
real(8), dimension(kms:kme+1), intent(in)  :: z
real(8), intent(in)                        :: hl, wl, rl
real(8), dimension(kms:kme), intent(out)   :: canop_a, canop_m
!- local variables
real(8), dimension(kms:kme)     :: pw
integer ::  k

do k = kms, kme
   if(k.le.kce-1) pw(k) = 1.d0
   if(k.ge.kce) pw(k) = 0.d0
end do

if(sf_surface_physics==1) then

   canop_a(:) = 0.0d0
   canop_m(:) = 1.0d0

else
   do k = kms, kme

      if((wl.eq.0.d0).or.(rl.eq.0.d0)) then
	       canop_a(k) = 0.d0
         canop_m(k) = 1.d0
      else
         canop_a(k) = wl*pw(k) / ((wl+rl)**2 - (wl**2)*pw(k))
         canop_m(k) = 1.d0 - (wl**2)*pw(k)/((wl+rl)**2)
      endif

   end do
end if

return
end subroutine canop
!=============================================================









!==========================================================
 subroutine ts_all(kce, tssl, tgsl, ts_r, ts_g, tw_r, tw_g, ts, tw )

implicit none
integer, intent(in)                             :: kce
real(8), dimension(4,kce-2), intent(out)        :: ts
real(8), dimension(4,kce-2, BLNU), intent(out)  :: tw
real(8), intent(out)                        :: tssl, ts_r, ts_g
real(8), dimension(nw), intent(out)         :: tgsl
real(8), dimension(RLNU), intent(out)       :: tw_r
real(8), dimension(GLNU), intent(out)       :: tw_g

ts_g    = 0.d0
ts_r    = 0.d0
tssl    = 0.d0
tgsl(:) = 0.d0
tw_r(:) = 0.d0
tw_g(:) = 0.d0
ts(:,:) = 0.0d0
tw(:,:,:) = 0.0d0

return
end subroutine ts_all
!===========================================================









!==========================================;
subroutine assign_urb_parm ( UTYPE )
implicit none
integer, intent(IN)      :: UTYPE

ZR           = ZR_TBL                 (UTYPE)
ROOF_WIDTH   = ROOF_WIDTH_TBL (UTYPE)
ROAD_WITH    = ROAD_WITH_TBL(UTYPE)
AH           = AH_TBL(UTYPE)
FRC_URB      = FRC_URB_TBL(UTYPE)
CAPR         = CAPR_TBL(UTYPE)
CAPB         = CAPB_TBL(UTYPE)
CAPG         = CAPG_TBL(UTYPE)
AKSR         = AKSR_TBL(UTYPE)
AKSB         = AKSB_TBL(UTYPE)
AKSG         = AKSG_TBL(UTYPE)
ALBR         = ALBR_TBL(UTYPE)
ALBB         = ALBB_TBL(UTYPE)
ALBG         = ALBG_TBL(UTYPE)
Z0B          = Z0B_TBL(UTYPE)
Z0G          = Z0G_TBL(UTYPE)
Z0R          = Z0R_TBL(UTYPE)
TRLEND       = TRLEND_TBL(UTYPE)
TBLEND       = TBLEND_TBL(UTYPE)
TGLEND       = TGLEND_TBL(UTYPE)
DDZR         = DDZR_TBL(UTYPE)
DDZB         = DDZB_TBL(UTYPE)
DDZG         = DDZG_TBL(UTYPE)
AHDIUPRF(:)  = AHDIUPRF_TBL(:)


THERMAL_INSOL_ROOF  = THERMAL_INSOL_ROOF_TBL
THERMAL_INSOL_WALL  = THERMAL_INSOL_WALL_TBL
RLNU         = RLNU_TBL
BLNU         = BLNU_TBL
GLNU         = GLNU_TBL
BOUNDR       = BOUNDR_TBL
BOUNDB       = BOUNDB_TBL
BOUNDG       = BOUNDG_TBL
AHOPTION     = AHOPTION_TBL
FRCURBOPT    = FRCURBOPT_TBL


if(.TRUE.) then
   write(20,*) ""
   write(20,*) "  *** Urban parameter ***"
   write(20,'(a20, f7.2)') "ZR", ZR
   write(20,'(a20, f7.2)') "ROOF_WIDTH", ROOF_WIDTH
   write(20,'(a20, f7.2)') "ROAD_WITH", ROAD_WITH
   write(20,'(a20, f7.2)') "AH", AH
   write(20,'(a20, f7.2)') "FRC_URB", FRC_URB
   write(20,'(a20, f10.0)') "CAPR", CAPR
   write(20,'(a20, f10.0)') "CAPB", CAPB
   write(20,'(a20, f10.0)') "CAPG", CAPG
   write(20,'(a20, f7.2)') "AKSR", AKSR
   write(20,'(a20, f7.2)') "AKSB", AKSB
   write(20,'(a20, f7.2)') "AKSG", AKSG
   write(20,'(a20, f7.2)') "ALBR", ALBR
   write(20,'(a20, f7.2)') "ALBB", ALBB
   write(20,'(a20, f7.2)') "ALBG", ALBG
   write(20,'(a20, f7.2)') "Z0B", Z0B
   write(20,'(a20, f7.2)') "Z0G", Z0G
   write(20,'(a20, f7.2)') "Z0R", Z0R
   write(20,'(a20, f7.2)') "TRLEND", TRLEND
   write(20,'(a20, f7.2)') "TBLEND", TBLEND
   write(20,'(a20, f7.2)') "TGLEND", TGLEND
   write(20,'(a20, f7.2)') "DDZR", DDZR
   write(20,'(a20, f7.2)') "DDZB", DDZB
   write(20,'(a20, f7.2)') "DDZG", DDZG
   write(20,'(a20,24f4.1)') "AHDIUPRF", AHDIUPRF
   write(20,'(a20, i4)') "THERMAL_INSOL_ROOF", THERMAL_INSOL_ROOF
   write(20,'(a20, i4)') "THERMAL_INSOL_WALL", THERMAL_INSOL_WALL
   write(20,'(a20, i4)') "RLNU", RLNU
   write(20,'(a20, i4)') "BLNU", BLNU
   write(20,'(a20, i4)') "GLNU", GLNU
   write(20,'(a20, i4)') "BOUNDR", BOUNDR
   write(20,'(a20, i4)') "BOUNDB", BOUNDB
   write(20,'(a20, i4)') "BOUNDG", BOUNDG
   write(20,'(a20, i4)') "AHOPTION", AHOPTION
   write(20,'(a20, i4)') "FRCURBOPT", FRCURBOPT
   write(20,*) ""
end if

end subroutine assign_urb_parm














!================================================================
subroutine read_landuse (luind)

implicit none
integer, intent(in)    :: luind
integer  :: IOSTATUS, allocate_status
integer  :: nuvar
character(len=300)  :: string
character(len=50) :: name
integer :: indx
character(LEN=100) :: a1
integer  :: i1
real(8)  :: x1, x2, x3, x4, x5, x6, x7

open ( 100, file='LANDUSE.TBL', access='SEQUENTIAL',        &
       status='OLD', action='READ', position='REWIND', &
       iostat=IOSTATUS )

if( IOSTATUS > 0 ) then
   stop
   print*, '   ERROR: open LANDUSE.TBL'
end if


nuvar = 0
READLOOP: do

   read(100,'(A300)', iostat=iostatus) string
   if (iostatus /= 0) exit READLOOP
   if (string(1:1) == "#") cycle READLOOP
   if (trim(string) == "") cycle READLOOP

   read(string,*,iostat=iostatus) i1, x1, x2, x3, x4, x5, x6, x7, a1

   if(i1==luind) then
      lu_slb%lu_id    = i1
      lu_slb%albd     = x1
      lu_slb%moist    = x2
      lu_slb%emiss    = x3
      lu_slb%z0       = x4
      lu_slb%bo       = x5
      lu_slb%lamda    = x6
      lu_slb%cs       = x7
      lu_slb%luname   = a1
      nuvar = 1
      print*, 'slab lu is: ',trim(lu_slb%luname)
   endif


enddo READLOOP
if(nuvar/=1) stop '  Error: lu index is not matching with LANDUSE.TBL: module_ini'

close(100)

end subroutine read_landuse

end module module_initialize
