!===============================
module module_sf_slab

use module_params, only: &
    history_interval,    &
    ganma,               &   ! lapse rate
    omega,               &   !
    sigma,               &   !
    g,                   &   ! gravity acceleration
    cp,                  &   !
    rd,                  &   !
    lv,                  &   ! latent heat
    tb,                  &   !
    pb,                  &   !
    dt,                  &   !
    kms, kme,            &   !
    soil_temp,           &   !
    nw,                  &   !
    lu_slb,              &   !
    ncid

use module_vars, only :  &    !
    istep,               &    ! time step
    kce,                 &    ! number of canopy layers
    MOLS                      ! Monin-Obukhov length (slab)
use module_phys_funcs
use module_io

!============
contains

!=======================================================
subroutine sf_slab  (                              &
                      nu_slab,                     & !- in
                      Z_T, u, v,                   & !- in
                      sw, lw,                      & !- in
                      t, ts, tw_g,                 & !- inout
                      qv,                          & !- in
                      tau_u, tau_v, tau_t, tau_qv  & !- out
                     )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
integer, intent(in)                           :: nu_slab
real(8), dimension(kms:kme), intent(in)       :: Z_T
real(8), dimension(kms:kme), intent(in)       :: u, v
real(8), intent(in)                           :: sw, lw
real(8), intent(inout)                        :: ts
real(8), dimension(kms:kme), intent(in)       :: t, qv
real(8), dimension(nw), intent(inout)         :: tw_g
real(8), intent(out)        :: tau_u, tau_v, tau_t, tau_qv

!- local variables
real(8)                        :: u1, v1, sq_u, UA
real(8)                        :: t_rib, tc1, tc, p
real(8)                        :: albedo, bo, lamda, cs, z0
real(8)                        :: ts1, ts2, d1, d2
real(8)                        :: rho, rib, ZA
real(8)                        :: rnet, h, le, g0, rdown, ta, tg
real(8)                        :: ch, cm, tau
real(8)                        :: hc, gc, lec, beta
real(8)                        :: gan_t
real(8)                        :: q, r
real(8)                        :: a, d, e, x, qv_sat, qva
integer                        :: ground_scheme, k_ref
real(8), dimension(nw)         :: kw, dx
real(8), parameter :: eps = 0.98d0   ! rate of longwave radiation
real(8)    :: soil_depth = 2.d0   ! soil layer depth (m)
real(8)    :: XXX, ALPHA, B1
!-----------------------------------------------------------------------------
if(istep==1) print*, "SLAB INDEX  : ", NU_SLAB
if(istep==1) print*, "SLAB ALBEDO : ", LU_SLB%ALBD
if(istep==1) print*, "SLAB Z0     : ", LU_SLB%Z0
if(istep==1) print*, "SLAB CS     : ", LU_SLB%CS
if(istep==1) print*, "SLAB LAMDA  : ", LU_SLB%LAMDA
if(istep==1) print*, "SLAB LAMDA  : ", LU_SLB%BO
if(istep==1) print*, "SLAB INDEX: ", NU_SLAB

k_ref = KMS + 1
! the soil parameter is set up according to land type
! the thermal properties of soil include:
!  albedo, i.e. reflectance of soil surface
!  bo, i.e bowen ration of soil surface
!  z0, i.e roughness length of soil surface
!  cs, i.e heat capacity of soil surface
!  lamda, ie thermal conductivity
albedo =   lu_slb%albd * 0.01d0
bo     =   lu_slb%bo
z0     =   lu_slb%z0 * 0.01d0
cs     =   lu_slb%cs
lamda  =   lu_slb%lamda

if(istep == 1) print*,'----', z0, albedo
! slab model considers the soil as bulk layer
! with depth soil_depth; soil layer is divided into 12
! sub-layers with the same thickness.

dx(:)  =  soil_depth/12.d0   ![m]
kw(:)  =   lamda/cs  ! kw is diffusion coefficient of soil
! it is determined by fraction of lamda
! to heat capacity

gan_t       = -0.01d0 + ganma     ! lapse rate of normal temp [K/m]
p           =  pb                 ! base pressure [hPa]
tc          =  tb - 273.15d0      ! base temperature (C)

ZA  =  Z_T(k_ref)
!- Chu y t
rho = 1.2d0 ! p * 100.d0 / rd /( t(k_ref)  +  tb  +  ZA*gan_t)
!print*, rho
!- calculate square velocity
u1 = u(k_ref)
v1 = v(k_ref)

sq_u = u1**2 + v1**2
if(sq_u.lt.1.d-5) sq_u = 1.d-5
if(sq_u.gt.1.d4)  sq_u = 1.d4

UA = dsqrt( SQ_U )

!- Exchange coefficents of heat (ch) and momentum (cm) flux is calculated
!  by Louis (1979) scheme;
!  To calculte ch & cm, bulk Richardson number and reference height is needed.
!               The reference height from which exchange processes take place.
!
!- To calculate bulk Richardson number, we need:
!        air temperature (of ref point):
!        height (of ref point):
!        surface temperature:
!        wind speed (squared):

t_rib = ZA * ganma + t(k_ref)       !- real temp at (k_ref) ---tb
!  The Bulk Richardson Number (BRN) is a dimensionless number relating
!  vertical stability and vertical shear(generally, stability divided by shear).
!  It represents the ratio of thermally produced turbulence and
!                             turbulence generated by vertical shear.
!  Practically, its value determines whether convection is free or forced.
!  High values indicate unstable and/or weakly sheared environments;
!  low values indicate weak instability and/or strong vertical shear.
!  In general, minus BRN indicates instability, and plus BRN indicates
!  stability.
!if (istep == 171005) print*, rib
rib = (g/tb) * ( t_rib-ts ) * za / sq_u !/ ZA

call check_value ( istep, 'check rib', rib )

if(rib.le.-15.d0) rib = -15.d0
!call louis79_scheme (  rib, ZA, z0, ch, cm  )
B1 = DLOG(Z0/(Z0*0.1D0))/0.4
CALL mos( B1, ZA, Z0, UA, T_RIB, TS,  & !-in
          RIB, CM, CH, MOLS )           !-inout

!---Radiation
!   snet is short-wave radiation absorbed by surface
!   lnet is long-wave radiation absorbed by surface
!   rnet is net radiation,including both short and long
!
!beta   = 0.0d0 !0.3d0
HC     =  RHO * CP * CH *UA
LEC    =  RHO * LV * beta * CH * UA
GC     =  LAMDA/DX(1)

RDOWN  =  sw  * ( 1.0d0 - albedo ) + eps * lw
TA     =  t(k_ref) + ZA * gan_t
TG     =  tw_g(2)
QVA    =  qv(k_ref)

!  use newton method to calculate the surface temperature (ST)
!  Among the heat fluxes and air temperature used to calculate ST
!  only radiation-related variable,rdown, was updated to
!  current time-step, ta, tg, hc are old variables from previous
!  step. Therefore, the ST calculated here is just
!  semi-updated variable.

call sf_temp_newton_bowen_ad    (                        &
                                rdown, eps, ta, tg,      & ! in
                                hc, gc, bo, tb,          & ! in
                                ts, rnet, h, le, g0      & ! inout
                                )

! then calculate soil temperature
! the deep boundary condition of soil is fixed as soil_temp
! it is based on assumption that at underground deep enough, the soil
! temperature does not change during a year.

tw_g(1)  = ts
tw_g(nw) = soil_temp - tc

!TW_G(NW) = TW_G( NW-1 )   ! no flux inside

call heat_1d(dx, kw,tw_g)


!  momentum fluxes is directly proportional to velocity strength
!    and exchange coefficient between air and surface;
!    have plus sign because momentum fluxes appear to move from
!    air to surface, i.e the momentum fluxes is likely to be
!    absorbed by surface due to friction process.
TAU   = CM * UA**2
TAU_U = TAU * U1/UA
TAU_V = TAU * V1/UA

!  However, heat fluxes (in form of temperature fluxes) is directly
!    proportional to velocity strength, ie the exchange will be strong
!    under windy day, that is more heat fluxes will be brough to air.
!  On the other hand, heat fluxes is directly proportional to the
!    difference between surface and air temperature, ie there larger
!    that difference is, the more heat fluxes will be lifted up.
!  Furthermore, the exchange coefficient between air and surface plays
!    important role
TAU_T  =  - H / RHO / CP 
TAU_QV =  - LE / RHO / LV     ! [kg/kg x m/s]


call put_var_nc(istep, "BRiNu", rib)
call put_var_nc(istep, "CM_SLB", cm)
call put_var_nc(istep, "CH_SLB", ch)
call put_var_nc(istep, "TAU_U_SLB", tau_u)
call put_var_nc(istep, "TAU_V_SLB", tau_v)
call put_var_nc(istep, "TAU_T_SLB", tau_t)
call put_var_nc(istep, "RNET_SLB", rnet)
call put_var_nc(istep, "H_SLB", h)
call put_var_nc(istep, "LE_SLB", le)
call put_var_nc(istep, "GR_SLB", g0)
call put_var_nc(istep, "TS_SLB", ts + tc)

return
end subroutine sf_slab
!===========================================================



!===============================
end module module_sf_slab
!===============================
