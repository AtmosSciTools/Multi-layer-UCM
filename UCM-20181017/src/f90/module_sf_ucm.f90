!===============================
module module_sf_ucm


use module_params, only :      &
         history_interval,     &
         g,                    &  ! gravity
         lv,                   &  ! latent heat
         ganma,                &  !
         lapse_rate,           &  !
         omega,                &
         sigma,                &
         cp,                   &
         rd,                   &
         tb,                   &
         pb,                   &
         dt,                   &
         kms,                  &
         kme,                  &
         ncid

use module_vars, only:    &
         kce,               &
         istep,             &
         hl, wl, rl,        &    ! canopy geometry (height, width, road)
         z,                 &    ! vertical grid height
         dz,                &    ! grid distance
         z_t,               &    ! temperature grid
         dz_t,              &    ! temperature grid distance
         svf_r,             &    ! sky-view factor of roof
         svf_g,             &    ! sky-view factor of ground
         svf_w,             &    ! sky-view factor of wall
         vf_gw,             &    ! view factor of ground to wall
         vf_wg,             &    ! view factor of wall to ground
         vf_wwp,            &    ! view factor of wall to parallel wall
         vf_wwv,            &    ! view factor of wall to ortho wall
         !read urbanparams.tbl
         zr,                &    ! building height
         roof_width,        &    ! roof width
         road_with,         &    ! road width
         ah,                &    ! anthropogenic heat flux
         frc_urb,           &    ! fraction of urban area
         capr,              &    ! heat capacity of roof
         capb,              &    ! heat capacity of buiding (walls)
         capg,              &    ! heat capacity of ground
         aksr,              &    !
         aksb,              &    !
         aksg,              &    ! of ground
         albr,              &    ! albedo of roof
         albb,              &    !        of building (walls)
         albg,              &    !        of ground
         z0b,               &    ! roughness length of walls
         z0g,               &    !                  of ground
         z0r,               &    !                  of roof
         trlend,            &    !
         tblend,            &    !
         tglend,            &    !
         ddzr,              &    !
         ddzb,              &    !
         ddzg,              &    !
         thermal_insol_roof,                &  !
         thermal_insol_wall,                &  !
         rlnu,                              &  !
         blnu,                              &  !
         glnu,                              &  !
         boundr,                            &  !
         boundb,                            &  !
         boundg,                            &  !
         ahoption,                          &  !
         frcurbopt,                         &  !
         ahdiuprf,                          &  !
         ! radiation-solar
         cosz, omg, declin,                 &
         molr, molb, molg, molc

use module_io
use module_phys_funcs
implicit none
real(8)  :: x1, x2, x3
real(8), parameter :: therinsul_lamda  = 0.03d0
real(8), parameter :: therinsul_cap    = 0.02d6
real(8), dimension(1:6)     ::    swnet_ucm, lwnet_ucm, lwdown_ucm, lwup_ucm

!============
contains

!=================================================================
subroutine ucm_swave_wall(                               &
                            elevation, alpha_s, alpha_c, & !- in
                            sd,                          & !- in
                            sdr                          & !- out
                            )
!=================================================================
! to calculate direct short wave radiation in particular direction
!    sd is short wave radiation on horizontal plane
!    sde, sdw, sds, sdn is for east, west, south, north direction.
!..................................................................
implicit none
real(8), intent(in)    :: elevation, alpha_s, alpha_c
real(8), intent(inout)    :: sd
real(8), dimension(4), intent(out)  :: sdr
real(8)                :: x, y, time, sdtop, x1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(elevation <= 0.d0)  then
   sdr(:) = 0.d0
   sd     = 0.d0
else
   y   =  dtan(elevation)   ; if ( y < 0.01d0 ) y = 0.01d0
   x   =  1.0d0 / y
   sdr(1) = -sd * x * alpha_s
   sdr(2) =  sd * x * alpha_s
   sdr(3) = -sd * x * alpha_c
   sdr(4) =  sd * x * alpha_c
end if
   where(sdr<0.d0)  sdr=0.d0
return
end subroutine ucm_swave_wall
!==================================================================





!======================================================
subroutine ucm_roof (                   &
                      sw, lw,           & !- in
                      u, v, p, t,       & !- in
                      ts, tw_r,         & !- inout
                      qv,               & !- inout
                      qa, q, qev,       & !- out
                      rnet, h, le, gr   & !- out
                     )
implicit none
real(8), dimension(kms:kme), intent(in)  :: u, v, t, p, qv
real(8), dimension(kms:kme), intent(in)       :: qa
real(8), intent(in)                           :: sw, lw
real(8), intent(inout)                        :: ts
real(8), dimension(rlnu), intent(inout)       :: tw_r
real(8), dimension(kms:kme), intent(inout)    :: q, qev
real(8), intent(out)                    :: h, rnet, gr, le

!- local variable
integer                        :: k_ref, k
real(8)                        :: gan_t, tc
real(8)                        :: albedo, bo
real(8)                        :: z0, rho, pa
real(8)                        :: ua, sq_u
real(8), dimension(rlnu)       :: kw, dx
real(8), parameter             :: eps = 0.9d0
real(8)                        :: ch, cd, rib, z_ref, za
real(8)                        :: rdown, ta, tg, qva
real(8)                        :: hc, gc, lec
real(8)                        :: beta
real(8)                        :: b1
!................................................................
! sensible heat is the product of heat-transfer coefficient and
! the difference between surface temp and reference air temp
! albr : albedo of roof surface
! bo: bowen ratio
! z0r: rough length of roof
if(istep==1) write(*,'(a,f7.1)') 'roof thickness', ddzr

dx(:)    =  ddzr/dfloat(rlnu-1)
albedo   =  albr
bo       =  5.d0
z0       =  z0r
kw(:)    =  aksr / capr
gc       =  aksr / dx(1)

! set insolation for roof
if(thermal_insol_roof==1) then
   if(istep==1) write(20,*) "-----insulation roof"
   do k=1,5
      kw(k) = therinsul_lamda / therinsul_cap
   end do
      gc    = therinsul_lamda / dx(1)
end if

gan_t    = lapse_rate + ganma   ![k/m]
k_ref    = kce
tc       = tb - 273.15d0

pa     =  p(k_ref)     ![hpa]
z_ref  =  z_t(k_ref)   !- height of referent point
ta     =  t(k_ref) + z_ref * gan_t

! air density (rho) = air pressure (pa) / rd / temp (kelvin)
rho    =  pa * 100.d0 / rd / (ta + tb )   !- air density
sq_u   = u(k_ref-1)**2 + v(k_ref-1)**2
if(sq_u.lt.1.d-5) sq_u = 1.d-5
ua     = dsqrt( sq_u )

za     =  dz(k_ref)/2.0d0     !- length from surface to referent point
rib    = (g/tb)*( ta - ts )* za / sq_u
!call louis79_scheme ( rib, za, z0, ch, cd  )
b1 = dlog(z0/(z0*0.1d0))/0.4
call mos( b1, za, z0, ua, ta, ts,  &           !-in
          rib, cd, ch, molr )                  !-inout

hc     =  rho * cp * ch * ua

!>>>> - newton method
rdown  =  sw * ( 1.0d0 - albedo ) + eps * lw
tg     =  tw_r(2)
qva    =  qv (k_ref)

call sf_temp_newton_bowen_ad    (                        &
                                rdown, eps, ta, tg,      & ! in
                                hc, gc, bo, tb,          & ! in
                                ts, rnet, h, le, gr      & ! inout
                                )

!>>>> step 5: calculate inside wall temperature
! set boundary condition

tw_r(1)  = ts
if( boundr == 1 ) then
   tw_r(rlnu) = tw_r(rlnu-1)
else
   tw_r(rlnu) = trlend  - tc
end if

call heat_1d_ucm( rlnu, dx, kw, tw_r)
!>>>> step four : calculate heat released by roof surface
! calculate sensible heat from roof of urban canopy
! h: w/m2 ; rho: kg/m3 ; cp: j/kgk
! q: k/s
q(k_ref) =  ((h+qa(kce)) /(rho*cp) )*( wl**2 / &
                (dz(k_ref)*(wl+rl)**2) )         ![ k/s ]
! specific humidity
! [j/kg x kg/m3 x kg/kg x m/s] = [w/m2]
qev(k_ref) = ((  wl/(wl+rl)  )**2) *le / (rho*lv)/dz(k_ref)
swnet_ucm(1)   = sw * ( 1.0d0 - albedo )
lwdown_ucm(1)  = eps * lw
lwup_ucm(1)    = eps * sigma*(ts+tb)**4
lwnet_ucm(1)   = lwdown_ucm(1) - lwup_ucm(1)
return
end subroutine ucm_roof
!===================================================================








!==================================================================
subroutine ucm_ground  (                                   &
                            u, v, p,                       & !- in
                            sw, lw, sd, ss, sdr,           & !- in
                            ratio_g, ratio,                & !- in
                            svf_g, vf_gw, hour,            & !- in
                            ts, t, ts_g, tw_g,             & !- inout
                            qv,                            & !- inout
                            tau_u, tau_v, tau_t, tau_qv,   & !- out
                            h, rnet, gr, le                & !- out
                          )
implicit none
real(8), dimension(kms:kme), intent(in)    :: u, v, p
real(8), intent(in)                        :: sw, lw, sd, ss
real(8), dimension(4), intent(in)          :: sdr ! short wave
real(8), intent(in)        :: ratio_g, svf_g
real(8), dimension(4, kce-2), intent(in) :: ratio     !  sunlit wall,
real(8), dimension(1:kce-2), intent(in)  :: vf_gw  ! vf ground to wall
integer, intent(in)                      :: hour
real(8), dimension(4, kce-2), intent(in) :: ts
real(8), dimension(kms:kme), intent(inout)   :: t, qv
real(8), dimension(glnu), intent(inout)      :: tw_g
real(8), intent(inout)                       :: ts_g
real(8), intent(out)           :: tau_u, tau_v, tau_t, tau_qv
real(8), intent(out)           :: h, rnet, gr, le
!- local variables
integer                        :: id, k
real(8)                        :: pa, rho, beta
real(8)                        :: u1, v1, sq_u, ua, za
real(8)                        :: albedo, bo, z0, tc, gan_t
real(8)                        :: rib
real(8)                        :: sdown, ldown  !, vfgw
real(8)                        :: lwall, swall, rdown
real(8)                        :: ch, cd
integer                        :: k_ref
real(8)                        :: hc, gc, lec
real(8)                        :: ta, tg, qva
real(8), dimension(glnu)       :: kw, dx           ! road diffusivity
real(8)                        :: b1, ahh
real(8), parameter             :: eps = 0.9d0
!--------------------------------------------------------------------

if(istep==1)  print*, "anthropogenic heat flux: ", ah

ahh = ah*ahdiuprf(hour+1)

k_ref   = kms + 1
if(istep==1)  write(*,'(a,f7.1)') 'ground thickness', ddzg

dx(:)       =  ddzg/dfloat(glnu-1)
gan_t       =  lapse_rate + ganma         ! [k/m]
pa          =  p(k_ref)                      ! [hpa]
tc          =  tb - 273.15d0      ! base temperature in c degree
albedo      =  albg
bo          =  3.d0
z0          =  z0g
ta          =  t(k_ref)
za          =  z_t(k_ref)

rho         = pa * 100.d0/ rd/ (ta + tb +  za * gan_t)
!- ground heat: kw is heat diffusivity
!= the product of thermal conductivity and capacity
kw(:)       =  aksg / capg

!- calculate square velocity
u1    =  u(k_ref)
v1 = v(k_ref)
sq_u  =  u1**2 + v1**2
if(sq_u.lt.1.d-5) sq_u = 1.d-5
ua    =  dsqrt( sq_u )
!- exchange coefficient
rib   =  (g/tb) * za * (ta - ts_g) / sq_u
b1    =  dlog(z0/(z0*0.1d0))/0.4
call mos( b1, za, z0, ua, ta, ts_g,  & !-in
              rib, cd, ch, molg )      !-inout

! >>>  : calc. heat budget at road surface
! reflected radiation from walls: loop id for 4 directions wall,
!                                 loop k for canopy layers
swall = 0.0d0
do id = 1, 4
  do k = 1, kce-2
     swall =   swall                                    &
            + albb  * vf_gw(k) * sdr(id) * ratio(id, k) &  ! direct
            + albb  * vf_gw(k) * ss                ! diffused
  end do
end do

!-total short wave equal sum of direct and diffused radiation on ground &
!                                          reflectd radiation from walls
sdown  = ( 1.0d0 - albedo ) * ( ratio_g * sd  +   svf_g * ss  +  swall )
!-total long wave released from walls to ground
lwall = 0.0d0
do id = 1, 4
   do k=1, kce-2
      lwall = lwall + vf_gw(k) * eps *  sigma  * ( tb + ts(id,k) ) **4
   end do
end do
!if(4.d0 * sum( vf_gw(:)) +svf_g .ne. 1.d0 ) then
!  print*, 4.d0 * sum(vf_gw(:))+svf_g
!  stop 'view factor ground'
!end if
ldown  = eps * ( svf_g * lw  +  lwall  )         ! absorbed by ground

hc     =  rho * cp * ch *ua       ! heat transfer coeffient ( w/m2k)
gc     =  aksg / dx(1)             ! heat transfer coefficient (w/km2)

! >>>     solve ts by newton method
rdown = sdown + ldown
tg    =  tw_g(2)
qva   =  qv(k_ref)

call sf_temp_newton_bowen_ad (                        &
                              rdown, eps, ta, tg,      & ! in
                              hc, gc, bo, tb,          & ! in
                              ts_g, rnet, h, le, gr    & ! inout
                              )
! >>>   solve wall temperature by 1d heat equation
!- bc of walls
tw_g(1)    =  ts_g
tw_g(glnu) =  tw_g(glnu-1)      !tglend  -  tc
call heat_1d_ucm ( glnu, dx, kw,tw_g)
tau_u      =     cd * ua  * u1               ! (m2/s2)
tau_v      =     cd * ua  * v1               ! cm * ua**2 is tau
tau_t      =   - h / rho / cp - ahh/rho/cp
tau_qv     =   - le / rho / lv                  ! [kg/kg x m/s]

swnet_ucm(2)   = sdown
lwdown_ucm(2)  = ldown
lwup_ucm(2)    = eps * sigma*(ts_g+tb)**4
lwnet_ucm(2)   = lwdown_ucm(1) - lwup_ucm(1)
return
end subroutine ucm_ground
!====================================================================






!=====================================================================
subroutine mucm_wall (                                     &
                       u, v, t,                            & !- in
                       sw, sd, ss, lw, sdr,                & !- in
                       ratio_g, ratio,                     & !- in
                       ts_g, ts, tw,                       & !- inout
                       qv,                                 & !- inout
                       qa, q, qev,                         & !- out
                       h, rnet, gr, le                     & !-out
                      )
implicit none
!- inout variables
real(8), dimension(kms:kme), intent(in)       :: u, v, t, qv
real(8), intent(in)                           :: sw, sd, ss, lw
real(8), intent(in)                           :: ts_g
real(8), dimension(4), intent(in)             :: sdr

real(8), dimension(4,kce-2), intent(in)           :: ratio
real(8), intent(in)                               :: ratio_g
real(8), dimension(kms:kme), intent(in)           :: qa
real(8), dimension(kms:kme), intent(inout)        :: q, qev
real(8), dimension(4,kce-2), intent(inout)        :: ts
real(8), dimension(4,kce-2,blnu), intent(inout)   :: tw
real(8), dimension(1:4,1:kce-2), intent(out)      :: h, rnet, gr, le
!- local variables
real(8), parameter             :: eps = 0.9d0
real(8), dimension(4)          :: tsw, xxx
real(8)                        :: u_cap, v_cap, q_cap
integer                        :: k , i, id, j
real(8)                        :: dzcap
real(8)                        :: gan_t, bo, p, c, albedo

real(8), dimension(4, kce-2)          :: hflux
real(8), dimension(4, kce-2)          :: snet, lnet, ldown, lup
real(8), dimension(4)                 :: sd_tag, sd_p, sd_v1, sd_v2
real(8), dimension(4, kce-2) :: ratio_tag, ratio_p, ratio_v1, ratio_v2
real(8), dimension(4, kce-2)          :: ts_tag, ts_p, ts_v1, ts_v2
real(8), dimension(4, kce-2, blnu)    :: tw_tag
real(8), dimension(blnu)              :: tw1d, kw, dx
real(8), dimension(kce-2)             :: z_cap
! - newton method
real(8)                               :: a, d, e, x, gc
! - for check
real(8)                               :: tc, tsg, large_u, hc, rho
real(8)                               :: s1, s2, s3, l1, l2, l3
real(8)                               :: x1, x2, x3, x4, x5

real(8)                               :: rdown, ta, tg, qva , xxx1
real :: xxx2
real(8) :: ahh
!------------------------------------------------------------------
ahh = 0.d0
if(istep==1) write(*,'(a,f7.1)') 'wall thickness', ddzb
dx(:)       =  ddzb/dfloat(blnu-1)
gan_t       =  lapse_rate + ganma     ! [k/m]
bo          =  1000.d0                  ! bowen rate [ ]
p           =  pb                     ! pressure [hpa]
albedo      =  albb              ! albedo [ ]
! calculate canopy layer height mean point
z_cap(1:kce-2)  = z_t(kms+1:kce-2+1)
tc              = tb - 273.15d00
kw(:)           = aksb / capb
gc              = aksb / dx(1)

if(thermal_insol_wall==1) then
if(istep==1) write(20,*) '---insulation walls'
   do k = 1,5
      kw(k) = therinsul_lamda / therinsul_cap
   end do
      gc    = therinsul_lamda / dx(1)
end if

!>>>>  - determine dimension variables  !okkkk
sd_tag(:)    = sdr(:)
sd_p(1)      = sdr(2)   ! sd_p is sd from parallel wall
sd_p(2)      = sdr(1)   ! 1 is east, 2 is west, 3 is south,
sd_p(3)      = sdr(4)   ! 4 is north
sd_p(4)      = sdr(3)
!- determine ratio under solar

ratio_tag        = ratio
ratio_p(1, :)    = ratio(2, :)
ratio_p(2, :)    = ratio(1, :)
ratio_p(3, :)    = ratio(4, :)
ratio_p(4, :)    = ratio(2, :)

!- determine  wall surface temp variables ( two dims )
!- then dim 1 is walls derection ; dim 2 is canopy layer nu

ts_tag      =  ts
ts_p(1, :)  =  ts(2, :)
ts_p(2, :)  =  ts(1, :)
ts_p(3, :)  =  ts(4, :)
ts_p(4, :)  =  ts(3, :)
tsg         =  tb + ts_g   ! (in kelvin) use old ts_ground

!- walls inside temperature
tw_tag(:,:,:)      = tw(:,:,:)
!---------------------------
!- loop from here
!---------------------------
do id = 1, 4            ! 1-east, 2-west, 3-south, 4-north
   do i = 1, kce-2      ! loop for canopy layer
      !>>>> step 2: calculate short wave radiation        [ w/m2 ]
      ! net short wave = short wave from
      !                  1) road 2) other wall and 3) solar
      ! s1 = reflected radiation: from road and wall (parallel)
      s1 = vf_wg(i)  * ratio_g * sd * albg   ! from road

      do k=1, kce-2
         s1  =  s1 + albedo *  (vf_wwp(i,k)+2.d0*vf_wwv(i,k)) * &
             ratio_p(id, k)  * sd_p(id)
      end do

      ! s2 is direct short way from atmosphere (only area exposed to sun)
      s2  =  ratio_tag(id, i) * sd_tag(id)
      ! s3 is diffusive short way from atmosphere
      s3  =  svf_w(i) * ss
      ! snet = ( 1-albedo ) absorb rate x total coming radiation
      snet(id, i) = ( 1.0d0 - albedo ) * ( s1 + s2 + s3 )

      !>>>> step 3: calculate long wave radiation                [ w/m2 ]

      l1 = eps * sigma * vf_wg(i) * tsg**4     ! from road
      do k=1, kce-2
         l1  =  l1 +  eps * sigma *  (vf_wwp(i, k)+2.d0*vf_wwv(i,k)) *   &
                                        ( tb + ts_p (id, k) ) **4
      end do

      l2 =  svf_w(i) * lw

      l3 =  eps * sigma * ( tb  +  ts_tag(id, i) ) ** 4

      ldown(id, i) = l1 + l2
      lup(id, i)   = l3
      lnet(id, i) =  ldown(id, i) - lup(id, i)
      rnet(id, i) = snet(id, i) + lnet(id, i)

      ! calculate sensible heat flux by jurjes formula        [ w/m2 ]

      large_u    = ( u(i+1) **2 + v(i+1) **2 ) ** 0.5d0
      if( large_u <= 5.0d0 ) then
         hc = 6.15d0 + 4.18d0 * large_u
      else
	       hc = 7.51d0 * large_u **0.78d0
      end if

      ! newton method
      a =  - eps * sigma
      d =  - (1.0d0 + 1.0d0/bo ) * hc   -  gc
      e =  snet(id, i) + ldown(id, i) +  gc * ( tw_tag(id,i, 2) + tb ) +  &
                (1.0d0 + 1.0d0/bo ) * hc * ( t(i+1)  + z(i+1) * ganma  +  tb  )

      !  call newton ( a, d, e, x )
      rdown  =  snet(id, i) + ldown(id, i)
      ta     =  t(i+1)  +  z(i+1) * gan_t
      tg     =  tw_tag(id,i, 2)
      qva    =  qv(i)

      call sf_temp_newton_bowen ( rdown, ta, tg, hc, gc, bo, x  )

      call sf_temp_newton_bowen_ad (                     &
                              rdown, eps, ta, tg,        & !-in
                              hc, gc, bo, tb,            & !-in
                              ts_tag(id,i), rnet(id,i),  & !-io
                              h(id,i),le(id,i),gr(id,i)  & ! inout
                              )
      !ts_tag(id, i) = x - tb
      ! inside_wall temperature

      do j=1, blnu   ; tw1d(j)  = tw_tag(id,i, j)   ;end do

      tw1d(1)  =   ts_tag(id, i)
      if( boundb == 1 ) then
          tw1d( blnu )=tw1d(blnu-1)
      else
          tw1d(blnu) = tblend - tc
      end if
      call heat_1d_ucm(blnu, dx, kw, tw1d)

      !- return values
      do j=1, blnu   ;tw_tag(id,i,j) = tw1d(j)    ;end do

      do j=1, blnu
          tw(id, i, j) = tw_tag(id,i, j)
      end do

      ts(id,i) =  ts_tag(id, i)
      call check_value ( istep, "ts wall", ts(id,i) )
      !>>>>>>>>>>>>>> calculate flux from new ts
      !------
      h(id, i)    =   hc *( ts_tag(id, i) - t(i+1) )       ! w/m2

      x1     =  t(i+1) + tb
      rho    =  p * 100.0d0 / rd / x1               ! equation of state
      hflux(id, i) =  h(id, i)/(cp*rho)

      le(id, i) = h(id, i)/bo                         !- latent heat [w/m2]
      gr(id, i) = gc * ( ts_tag(id, i) - tw_tag(id,i,2) )
   end do   ! end loop directions
end do

!...............................................
!   end loop
!..............................................

!- calculate heat from canyon
do i=1, kce-2
   x1 =  sum(hflux(:,i)) * dz(i+1) *wl              ![km3/s]
   x2 = x1 / ((rl*rl + 2.0d0 * rl*wl)*dz(i+1))      ![k/s]
   q(i+1) = x2
end do

! km/s
q(kms+1) = q(kms+1) + ahh/(cp*rho) / dz(kms+1)

x1=1.0d0/dfloat(kce-2)

swnet_ucm(3:6)   = sum(snet,dim = 2)   / dfloat(kce-2)
lwdown_ucm(3:6)  = sum(ldown,dim = 2)  / dfloat(kce-2)
lwup_ucm(3:6)    = sum(lup,dim = 2)    / dfloat(kce-2)
lwnet_ucm(3:6)   = lwdown_ucm(3:6) - lwup_ucm(3:6)


return
end subroutine mucm_wall
!================================================================










!============================================
subroutine heat_1d_ucm (nw, dx, kw, tw)
use module_params, only : dt

!  solve heat conduction eq. by 1 dimension
implicit none
integer, intent(in)                    :: nw
real(8), dimension(nw), intent(in)     :: dx
real(8), dimension(nw), intent(inout)  :: tw, kw
!- local vars
real(8),dimension(nw)         :: dxx, a, b, c, d, gg, ss
integer                                :: k


do k = 2, nw-1
   a(k) = - kw(k-1) / ( dx(k) * (dx(k) + dx(k+1) ) * 0.5d0 ) * dt
   c(k) = - kw(k)   / (dx(k+1)* (dx(k) + dx(k+1) ) * 0.5d0 ) * dt
   b(k)   =   1.d0 - a(k)  - c(k)
   d(k)   =   tw(k)
end do

d(2)      =  d(2)     -  a(2)*tw(1)
a(2)      =  0.d0
d(nw-1)   =  d(nw-1)  -  c(nw-1)*tw(nw)
c(nw-1)   =  0.d0

!- solve tri
!- initialize g and s
gg(2) = b(2)
ss(2) = d(2)
!- solve for vectors c-prime and d-prime
do k = 3,nw-1
   gg(k) = b(k) - a(k)*c(k-1)/gg(k-1)
   ss(k) = d(k) - a(k)*ss(k-1)/gg(k-1)
enddo
!- initialize x
   tw(nw-1) = ss(nw-1) / gg(nw-1)
!- solve for x from the vectors c-prime and d-prime
do k = nw-2, 2, -1
   tw(k) = ( ss(k) - c(k) * tw(k+1))/gg(k)
end do
return
end subroutine heat_1d_ucm
!===============================================================



!===============================
end module module_sf_ucm
!===============================
