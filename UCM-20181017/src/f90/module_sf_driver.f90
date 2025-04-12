module module_sf_driver

use module_params,    only : &
         kms, kme,           &
         nw,                 &
         sltype,             &
         lat, lon,           &
         sf_surface_physics, &
         sw_opt, lw_opt,     &
         sd_opt, ss_opt

use module_vars,    only :  &    !
         !istep,            &    ! time step
         kce,               &    ! number of canopy layers
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
         theta_t, qv_t,     &
         u_t, v_t,          &
         sw_b, lw_b,        &
         cpu_sf, cpu_rt,    &
         cputime

use module_sf_slab
use module_sf_ucm

contains
!=======================================================================
subroutine surface_driver (                                    &
                          hour, elevation, alpha_s, alpha_c,   & !- i
                          sw, sd, ss, lw,                      & !- i
                          u, v, p,                             & !- i
                          t, q, qa, qev, tssl, tgsl,           & !- io
                          ts_g, ts_r, tw_r, tw_g, ts, tw,      & !- io
                          qv, rho_a,                           & !- io
                          tau_u, tau_v, tau_t, tau_qv          & !- o
                          )
!=======================================================================
!- in case of urban canopy model
!  aim to calculate flux ( tau ) from surface and sensible heat ( q )
!  from surface of walls and roof of building.
implicit none
!- input variables
integer, intent(in)                         :: hour
real(8), intent(in)                         :: alpha_s, alpha_c
real(8), intent(in)                         :: elevation
real(8), intent(inout)                      :: sw, sd, ss, lw
real(8), dimension(kms:kme), intent(in)     :: u, v, p
!- inout variables
real(8), intent(out)                        :: tau_u, tau_v, tau_t
real(8), intent(out)                        :: tau_qv
real(8), dimension(kms:kme), intent(inout)  :: t, qv, rho_a
real(8), dimension(kms:kme), intent(inout)  :: qev, q, qa
real(8), intent(inout)                      :: ts_r, ts_g, tssl

real(8), dimension(rlnu), intent(inout)            :: tw_r
real(8), dimension(glnu), intent(inout)            :: tw_g
real(8), dimension(1:4,1:kce-2), intent(inout)     :: ts
real(8), dimension(1:4,1:kce-2,1:blnu), intent(inout)    :: tw
real(8), dimension(nw),intent(inout)               :: tgsl

! dimension: roof, ground, east, west, south, north
!- local vars
integer                          :: k, n, m, nu_slab
real(8)                          :: xx
real(8)                          :: ratio_g
real(8), dimension(4, kce-2)     :: ratio
real(8), dimension(4)            :: sdr
real(8)                          :: tau_uu, tau_vv, tau_tt, tau_qq
real(8)                          :: alat, pi,elevation_s
real(8)               :: frac_urb
! for inout
real(8), dimension(1:6) :: rnet_ucm, h_ucm, le_ucm, gr_ucm, sd_ucm,    &
                           ts_ucm, cm_ucm, ch_ucm, ri_ucm, shade_ucm,  &
                           swnet, lwnet, lwdown, lwup
real(8), dimension(1:4,1:kce-2)    :: rnetw, hw, lew, grw
real(8)                 :: rnetg, rnetr, hr, hg, leg, ler, grr, grg
real(8) :: t1_rt, t2_rt, frc_bld
!-----------------------------------------------------------------------


pi = 4.d0*datan(1.d0)
alat      =  pi/180.d0 * lat
frc_bld   =  wl**2 / (wl+rl)**2

frac_urb  = (frc_urb - frc_bld) / ( 1.d0 - frc_bld )
nu_slab   = sltype
! step 1
! downward short wave solar radiation is divided into
! direct radiation (sd) and diffusive radiation (ss)
if (sw_opt > 0) then
   sw  =  sw_b
   sd  =  sw*0.8d0
   ss  =  sw - sd
endif
if (lw_opt > 0) lw  =  lw_b

call put_var_nc(istep, "SOL_ELEV", elevation)
call put_var_nc(istep, "AZTH_C", alpha_c)
call put_var_nc(istep, "AZTH_S", alpha_s)
call put_var_nc(istep, "SW", sw)
call put_var_nc(istep, "SD", sd)
call put_var_nc(istep, "SS", ss)
call put_var_nc(istep, "LW", lw)

surface_scheme: select case(sf_surface_physics)
case(1)                                          ! slab model
call sf_slab  (                               &
                 nu_slab,                     & !- in
                 z_t, u, v,                   & !- in
                 sw,lw,                       & !- in
                 t, ts_g, tgsl,               & !- inout
                 qv,                          & !- in
                 tau_u, tau_v, tau_t, tau_qv  & !- out
               )

case(2)
call sf_slab  (                               &
                 nu_slab,                     & !- in
                 z_t, u, v,                   & !- in
                 sw,lw,                       & !- in
                 t, tssl, tgsl,               & !- inout
                 qv,                          & !- in
                 tau_uu, tau_vv, tau_tt, tau_qq  & !- out
               )

!- step 1: calculate ratio under solar of each surface at each layer
!  by use building, road feature and solar condition
!  output: ration of each surface at each layer under solar
rnet_ucm  = 0.d0
h_ucm     = 0.d0
le_ucm    = 0.d0
gr_ucm    = 0.d0
ts_ucm    = 0.d0
cm_ucm    = 0.d0
ch_ucm    = 0.d0
ri_ucm    = 0.d0
shade_ucm = 0.d0
xx = ts_g
! to advoid using update ts_g to calculate walls surface temp. (ts).
! set xx as pre-update ts_g, then calculate new xx
call cpu_time(t1_rt)
call shade_multi (                                   &
                   elevation, alpha_s, alpha_c,      & !- in
                   ratio_g, ratio                    & !- out
                 )
call cpu_time(t2_rt)
cpu_rt = cpu_rt + t2_rt - t1_rt
!- step 2: calculate direct radiation at each direction
call ucm_swave_wall(                             &
                    elevation, alpha_s, alpha_c, & !- in
                    sd,                          & !- in
                    sdr                          & !- out
                    )

sd_ucm(1:2)  = sd
sd_ucm(3:6)  = sdr
!   call chck_rad ( istep, hl, wl, rl, sd, sdr, ratio, ratio_g )
!- step 3: calculate all surface temperature
!          calculate tau from road to first air layer
!          calculate heat flux q from each urban canopy layer.
!          calculate walls inside temperature. tw
call ucm_ground  (                                 &
                    u, v, p,                       & !- in
                    sw, lw, sd, ss, sdr,           & !- in
                    ratio_g, ratio,                & ! -in
                    svf_g, vf_gw, hour,            & !- in
                    ts, t, xx, tw_g,               & !- inout
                    qv,                            & !- inout
                    tau_u, tau_v, tau_t, tau_qv,   & !- out
                    hg, rnetg, grg, leg            & !- out
                  )
tau_t  = tau_t  * frac_urb   + tau_tt  * (1.d0 - frac_urb)
tau_u  = tau_u  * frac_urb   + tau_uu  * (1.d0 - frac_urb)
tau_v  = tau_v  * frac_urb   + tau_vv  * (1.d0 - frac_urb)
tau_qv = tau_qv * frac_urb   + tau_qq  * (1.d0 - frac_urb)
! update xx, but don't return xx to ts_g at this step.
call mucm_wall (                            &
                    u, v, t,                & !- in
                    sw, sd, ss, lw, sdr,    & !- in
                    ratio_g, ratio,         & !- in
                    ts_g, ts, tw,           & !- inout
                    qv,                     & !- inout
                    qa, q, qev,             & !- out
                    hw, rnetw, grw, lew     &
                 )
!  >>  all update surface temprerature
!- step 6: calculate sensible heat, exhausted from roof of building
call  ucm_roof(                           &
                     sw, lw,              & !- in
                     u, v, p, t,          & !- in
                     ts_r, tw_r,          & !- inout
                     qv,                  & !- inout
                     qa, q, qev,          & !- out
                     rnetr, hr, ler, grr )   !- out
q = q * frac_urb ! very important
ts_g = xx   ! finished calc. all walls sfctemp. return xx to ts_g

shade_ucm(2)    = ratio_g
shade_ucm(3:6)  = sum(ratio,dim=2) / dfloat(kce-2)

rnet_ucm(1)    = rnetr
rnet_ucm(2)    = rnetg
rnet_ucm(3:6)  = sum(rnetw,dim=2) / dfloat(kce-2)

h_ucm(1)    = hr
h_ucm(2)    = hg
h_ucm(3:6)  = sum(hw,dim=2) / dfloat(kce-2)

le_ucm(1)    = ler
le_ucm(2)    = leg
le_ucm(3:6)  = sum(lew,dim=2) / dfloat(kce-2)

gr_ucm(1)    = grr
gr_ucm(2)    = grg
gr_ucm(3:6)  = sum(grw,dim=2) / dfloat(kce-2)

ts_ucm(1)   = ts_r + tb - 273.15d0
ts_ucm(2)   = ts_g + tb - 273.15d0
ts_ucm(3:6) = sum(ts,dim = 2)  / dfloat(kce-2) + tb - 273.15d0


CALL put_var2d_nc ( istep, "RNET_UCM",  RNET_UCM  )
CALL put_var2d_nc ( istep, "H_UCM",     H_UCM     )
CALL put_var2d_nc ( istep, "LE_UCM",    CH_UCM    )
CALL put_var2d_nc ( istep, "GR_UCM",    GR_UCM    )
CALL put_var2d_nc ( istep, "TS_UCM",    TS_UCM    )
CALL put_var2d_nc ( istep, "SD_UCM",    SD_UCM    )
CALL put_var2d_nc ( istep, "CM_UCM",    CM_UCM    )
CALL put_var2d_nc ( istep, "CH_UCM",    CH_UCM    )
CALL put_var2d_nc ( istep, "SHADE_UCM", SHADE_UCM )

CALL put_var2d_nc ( istep, "SWNET",     SWNET_UCM    )
CALL put_var2d_nc ( istep, "LWNET",     LWNET_UCM    )
CALL put_var2d_nc ( istep, "LWDOWN",    LWDOWN_UCM   )
CALL put_var2d_nc ( istep, "LWUP",      LWUP_UCM     )

end select surface_scheme


return
end subroutine surface_driver
!=================================================================









!=================================================================
subroutine chck_rad ( istep, hl, wl, rl, sd, sdr, ratio, ratio_g )

implicit none
integer, intent(in)                      :: istep
real(8), dimension(4), intent(in)        :: sdr
real(8), dimension(4, kce-2), intent(in) :: ratio
real(8), intent(in)                      :: hl, wl, rl
real(8), intent(in)                      :: sd, ratio_g

integer                                  :: id
real(8), dimension(4)                    :: sum_rate
real(8)        :: sum_rad, wall_rate, road_rate, roof_rate

sum_rate = sum(ratio, dim=2)/dfloat(kce-2)
wall_rate = wl*hl/(wl+rl)**2
roof_rate = wl**2/(wl+rl)**2
road_rate = 1.0d0 - roof_rate

sum_rad = 0.d0
do id = 1, 4
   sum_rad = sum_rad + sdr(id) * sum_rate(id) * wall_rate
end do

sum_rad =  sum_rad + (roof_rate + ratio_g * road_rate)* sd

return
end subroutine chck_rad
!==================================











!==============================================================
subroutine shade_multi(                                  &
                         elevation, alpha_s, alpha_c,    & !-i
                         ratio_g, ratio                  & !-o
                       )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  calculate ratio of sunlit area by "ray tracing" method
!  output variables
!     ratio_g (scalar): for road
!     ratio (direction: e w s n, canop_layer: kce - 2)
!  input variables
!     elevation: solar elevation
!     alpha_s:   solar azimuth angle (sin)
!     alpha_c:   solar azimuth angle (cosin)
!     hl:        canopy height
!     wl:        building width
!     rl:        road width
!     z:         grid-point height
!     dz:        grid-point distance
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
real(8), intent(in)                 :: elevation, alpha_s, alpha_c
real(8), dimension(4,kce-2), intent(out)  :: ratio
real(8)                                   :: ratio_g
real(8), dimension(kce-2)   :: ratio_e, ratio_w, ratio_s, ratio_n
!- local variables
integer                      :: k
real(8)                      :: x, y, xa, ya
real(8)                      :: xs, ys
real(8)                    :: sq_shade, sq_road, sq_wall, sq_layer
real(8)                      :: shade_g
real(8)                      :: a0
real(8)                      :: dzcap, h1
real(8), dimension(kce-2)    :: z_cap
real(8), dimension(kce-2)    :: ash_e, ash_w, ash_s, ash_n
real(8), dimension(kce-2)    :: e1, e2, w1, w2, s1, s2, n1, n2
real(8), dimension(kce-2)    :: shade_p, shade_v
real(8)                      :: dsh1, dsh2, sum_sh
real(8)                           :: dsh, x1, x2
real(8), dimension(4)             :: rate
real(8)                           :: tanel


dzcap = hl/dfloat(kce-2)

z_cap(1) = 0.0d0
do k = 2, kce-2
   z_cap(k) = z_cap(k-1) + dzcap
end do

sq_layer = dzcap * wl
tanel    = dtan(elevation)  ; if(tanel<=0.1d0) tanel = 0.1d0

!- shadow area
x   =  hl*alpha_s/tanel
y   =  -hl*alpha_c/tanel

xa = dabs(x)
ya = dabs(y)

a0 = ( xa + ya ) * wl
!-before sun rise
if ( elevation.lt.0.0d0 ) then
   ratio(:,:) = 0.d0
   ratio_e(:) = 0.0d0
   ratio_w(:) = 0.0d0
   ratio_s(:) = 0.0d0
   ratio_n(:) = 0.0d0
   ratio_g    = 0.0d0
else
!- when solar is in south-east direction
!- shade is in north-west direction
   if( (x < 0.0d0).and. (y >= 0.0d0) ) then
      xs = xa
      ys = ya
      call shade_q(                              &
	            elevation, xs, ys,           & !- in
	            hl, wl ,rl,                  & !- in
	            sum_sh, shade_p, shade_v     & !- out
                   )

      dsh1 = sum_sh
      e1(:) = shade_p(:)
      s1(:) = shade_v(:)

      xs = ya
      ys = xa

      call shade_q (                              &
	             elevation, xs, ys,           & !- in
	             hl, wl ,rl,                  & !- in
	             sum_sh, shade_p, shade_v     & !- out
	            )

      dsh2 = sum_sh
      s2(:) = shade_p(:)
      e2(:) = shade_v(:)

      ash_e(:) = e1(:) + e2(:)
      ash_w(:) = sq_layer
      ash_s(:) = s1(:) + s2(:)
      ash_n(:) = sq_layer

      dsh = dsh1 + dsh2
   end if

   !- when solar is in north-east direction
   !- shade is in south-west direction

   if( (x < 0.0d0).and. (y < 0.0d0) ) then
      xs = xa
      ys = ya
      call shade_q(                              &
	            elevation, xs, ys,           & !- in
	            hl, wl ,rl,                  & !- in
	            sum_sh, shade_p, shade_v     & !- out
                   )

      dsh1 = sum_sh
      e1(:) = shade_p(:)
      n1(:) = shade_v(:)

      xs = ya
      ys = xa

      call shade_q(                              &
	            elevation, xs, ys,           & !- in
	            hl, wl ,rl,                  & !- in
	            sum_sh, shade_p, shade_v     & !- out
	           )

      dsh2 = sum_sh
      n2(:) = shade_p(:)
      e2(:) = shade_v(:)

      ash_e(:) = e1(:) + e2(:)
      ash_w(:) = sq_layer
      ash_s(:) = sq_layer
      ash_n(:) = n1(:) + n2(:)

      dsh = dsh1 + dsh2
   end if

   !- when solar is in north-west direction
   !- shade is in south-east direction
   if( (x >= 0.0d0).and. (y < 0.0d0) ) then
      xs = xa
      ys = ya
      call shade_q(                              &
	            elevation, xs, ys,           & !- in
	            hl, wl ,rl,                  & !- in
	            sum_sh, shade_p, shade_v     & !- out
                   )
      dsh1 = sum_sh
      w1(:) = shade_p(:)
      n1(:) = shade_v(:)

      xs = ya
      ys = xa

      call shade_q(                              &
	            elevation, xs, ys,           & !- in
                    hl, wl ,rl,                  & !- in
	            sum_sh, shade_p, shade_v     & !- out
	           )

      dsh2 = sum_sh
      n2(:) = shade_p(:)
      w2(:) = shade_v(:)
      ash_e(:) = sq_layer
      ash_w(:) = w1(:) + w2(:)
      ash_s(:) = sq_layer
      ash_n(:) = n1(:) + n2(:)

      dsh = dsh1 + dsh2


   end if

   !- when solar is in south-west direction
   !- shade is in north-east direction

   if( (x >=0.0d0).and. (y >= 0.0d0) ) then
      xs = xa
      ys = ya
      call shade_q(                              &
	            elevation, xs, ys,           & !- in
	            hl, wl ,rl,                  & !- in
	            sum_sh, shade_p, shade_v     & !- out
	           )

      dsh1 = sum_sh
      w1(:) = shade_p(:)
      s1(:) = shade_v(:)


      xs = ya
      ys = xa
      call shade_q(                              &
	            elevation, xs, ys,           & !- in
                    hl, wl ,rl,                  & !- in
	            sum_sh, shade_p, shade_v     & !- out
	           )

      dsh2 = sum_sh
      s2(:) = shade_p(:)
      w2(:) = shade_v(:)

      ash_e(:) = sq_layer
      ash_w(:) = w1(:) + w2(:)
      ash_s(:) = s1(:) + s2(:)
      ash_n(:) = sq_layer

      dsh = dsh1 + dsh2

   end if

   do k = 1, kce-2
     if(ash_e(k)<=0.d0) ash_e(k) = 0.d0
     if(ash_w(k)<=0.d0) ash_w(k) = 0.d0
     if(ash_s(k)<=0.d0) ash_s(k) = 0.d0
     if(ash_n(k)<=0.d0) ash_n(k) = 0.d0

     if(ash_e(k)>=sq_layer) ash_e(k) =sq_layer
     if(ash_w(k)>=sq_layer) ash_w(k) =sq_layer
     if(ash_s(k)>=sq_layer) ash_s(k) =sq_layer
     if(ash_n(k)>=sq_layer) ash_n(k) =sq_layer

     ratio_e(k) = 1.d0 - ash_e(k) / sq_layer  ; ratio(1, k) = ratio_e(k)
     ratio_w(k) = 1.d0 - ash_w(k) / sq_layer  ; ratio(2, k) = ratio_w(k)
     ratio_s(k) = 1.d0 - ash_s(k) / sq_layer  ; ratio(3, k) = ratio_s(k)
     ratio_n(k) = 1.d0 - ash_n(k) / sq_layer  ; ratio(4, k) = ratio_n(k)
   end do

   sq_road = 2.0d0*rl*wl + rl**2
   sq_shade = wl * ( dabs(xs) + dabs(ys) )

   shade_g = sq_shade - dsh
   if( shade_g >= sq_road ) shade_g = sq_road
   ratio_g  = 1.d0 - shade_g/sq_road
end if

rate = sum(ratio, dim=2) / dfloat(kce-2)

return
end subroutine shade_multi
!==========================================================








!==============================================================
 subroutine shade_q (                                  &
                          elevation, xs, ys,           & !- in
                          hl, wl ,rl,                  & !- in
                          sum_sh, shade_p, shade_v     & !- out
                     )
! xs, ys is shade length in x and y direction, respectly
! sum_sh is sum of shade area on the walls
! shade_p is
! shade_v is
implicit none
real(8), intent(in)            :: elevation, xs, ys, hl, wl, rl
real(8), intent(out)                   :: sum_sh
real(8), dimension(kce-2), intent(out) :: shade_p, shade_v
!- local
! bld is number of buildings
integer, parameter           :: wall_n = 10
integer, parameter           :: bld    = 5   ! reference building
integer                      :: i, j, k,  n, lgc1, lgc2
real(8)                      :: alpha_t
real(8)                      :: dw, dw1, dw2
real(8)                   :: x1, y1, z1, x2, y2, z2, xy, dz1, dz2
real(8)                      :: a0
real(8)                      :: hs_p, hs_v
integer, dimension(wall_n)   :: pv
real(8), dimension(wall_n)   :: wx, wy, a1, b1, c1
real(8), dimension(wall_n)   :: zsh, xsh, ysh, dsh, hsh
real(8), dimension(bld, 2)   :: a, b
real(8), dimension(kce-2)    :: hsh_p, hsh_v, z_cap, hs, s1, s2
real(8)                      :: dzcap, h1


alpha_t = ys / xs
a0 = xs * ys                     ! area of shade on the ground
lgc1 = 1                         ! logic coefficient 1
lgc2 = 1                         ! logic coefficient
xy = (xs**2 + ys**2)**0.5d0      !

dzcap = hl / dfloat(kce-2)    ! height of one canopy layer
z_cap(1) = 0.0d0              ! first z_cap height

do k = 2, kce-2
   z_cap(k) = z_cap(k-1) + dzcap     ! z_cap height
end do

dw    = wl /  dfloat(wall_n)     ! divide wl
dw1   = dw                       !
dw2   = dw / alpha_t             !
wy(1) = dw  /2.0d0               !

do i= 2, wall_n
   wy(i) = wy(i-1) + dw
   wx(i) = 0.0d0
end do
do i = 1, bld      ! ref. building number,1 parallel,2vertical
   a(i, 1) = rl + dfloat(i-1) * (wl+rl)
   a(i, 2) = a(i, 1) + wl
   b(i, 1) = dfloat(i-1) * (wl+rl)
   b(i, 2) = b(i, 1) + wl
end do

sum_sh = 0.0d0

do n = 1, wall_n                   ! start loop
   ! calculate parallel
    do i = 1, bld                 ! loop reference building number
       x1 = a(i, 1)
       y1 = wy(n) + x1 *alpha_t
       if(x1 > xs) then
          lgc1 = 0
          exit
       end if
       do j = 1, bld
	  if( (y1 >= b(j,1)).and.(y1 <= b(j,2))) then
	     lgc1  = 1
             exit
	  else
	     lgc1 = 0
	  end if
       end do

       if(lgc1==1) then
          z1 = (x1**2 + (x1*alpha_t)**2)**0.5d0
          dz1 = ( xs - x1 ) * dw
          exit
       end if

   end do

   ! calculate vertical
   do i = 1, bld
      y2 = b(i,1)
      x2 = ( y2 - wy(n) ) / alpha_t
      if(y2 > ys) then
         lgc2 = 0
         exit
      end if
      do j = 1, bld
	 if( (x2 >= a(j,1)).and.(x2 <= a(j,2))) then
	    lgc2  = 1
	    exit
	 else
	    lgc2 = 0
	 end if
      end do
      if(lgc2==1) then
         z2 = (x2**2 + (x2*alpha_t)**2)**0.5d0
	 dz2 = ( xs - x2 ) * dw
         exit
      end if
   end do

   if((lgc1==1).and.(lgc2==1)) then
      if(z1 < z2) then
	 zsh(n) = z1
         xsh(n) = x1
	 ysh(n) = y1
         dsh(n) = dz1
         pv(n)  = 1
      end if
      if(z1 >= z2) then
	 zsh(n) = z2
	 xsh(n) = x2
	 ysh(n) = y2
         dsh(n) = dz2
         pv(n)  = 2
      end if
   end if

   if((lgc1==0).and.(lgc2==1)) then
      zsh(n) = z2
      xsh(n) = x2
      ysh(n) = y2
      dsh(n) = dz2
      pv(n)  = 2
   end if

   if((lgc1==1).and.(lgc2==0)) then
      zsh(n) = z1
      xsh(n) = x1
      ysh(n) = y1
      dsh(n) = dz1
      pv(n)  = 1
   end if

   if((lgc1==0).and.(lgc2==0)) then
      zsh(n) = 0.0d0
      xsh(n) = 0.0d0
      ysh(n) = 0.0d0
      dsh(n) = 0.0d0
      pv(n)  = 0
   end if
   sum_sh = sum_sh + dsh(n)
end do



do n = 1, wall_n
   hsh(n) = hl -  zsh(n) * dtan(elevation)
end do
!- parallel walls and vertical walls
!- 1)  shade area by layer
shade_p(:) = 0.0d0
shade_v(:) = 0.0d0
do n = 1, wall_n
   do k = 1, kce-2
      h1 = hsh(n) - z_cap(k)
      if( h1 <= 0.0d0 ) hs(k) = 0.0d0
      if( (h1 > 0.0d0) .and. (h1 < dzcap) ) hs(k) = h1
      if( h1 >= dzcap ) hs(k) = dzcap
   end do
   !- parallel
   if(pv(n) == 1) then
      do k = 1, kce-2
         hsh_p(k) = hs(k)
         s1(k) = hsh_p(k) * dw1
         shade_p(k) = shade_p(k) + s1(k)
      end do
   end if
   !- vertical
   if(pv(n) == 2) then
      do k = 1, kce-2
         hsh_v(k) = hs(k)
         s2(k) = hsh_v(k) * dw2
         shade_v(k) = shade_v(k) + s2(k)
      end do
   end if

end do
return
end subroutine
!=========================




























!============================================================================
  subroutine shade(                                             &
                    istep,                                      & !- in
                    hl, wl, rl,                                 & !- in
                    elevation, alpha_s, alpha_c,                & !- in
                    ratio                                       & !- out
                  )
  implicit none

   integer, intent(in)                       :: istep
   real(8), intent(in)                       :: hl, wl, rl
   real(8), intent(in)                       :: elevation, alpha_s, alpha_c
   real(8), dimension(5), intent(out)         :: ratio
   real(8)        :: ratio_r, ratio_e, ratio_w, ratio_s, ratio_n

!- local variables

   real(8)    :: x, y
   real(8)    :: xs, ys
   real(8)    :: sq_shade, sq_road, sq_wall
   real(8)    :: hs_e, hs_w, hs_s, hs_n
   real(8)    :: ws_e, ws_w, ws_s, ws_n
   real(8)    :: shade_e, shade_w, shade_s, shade_n, shade_r
   real(8)    :: a0, a_shade


   if((hl.eq.0.d0).or.(wl.eq.0.d0).or.(rl.eq.0.d0)) then


        ratio_r = 1.d0
        ratio_e = 1.d0
        ratio_w = 1.d0
        ratio_s = 1.d0
        ratio_n = 1.d0
        x = 0.d0
        y = 0.d0

       else

        x =  hl*alpha_s/dtan(elevation)
        y =  -hl*alpha_c/dtan(elevation)


        a0 = ( x + y ) * hl

!-neu mat troi chua moc
      if(elevation.lt.0.d0) then

          hs_e = hl
          hs_w = hl
          hs_s = hl
          hs_n = hl

         else

!- - - - the hight of ( east , west )'s shade - - -

!- neu mat troi moc


         if((x.ge.0.d0).and.(x.lt.rl)) then
             hs_w = 0.d0
             hs_e = hl
         end if

         if(x.ge.rl) then
             hs_w = hl - rl*dtan(elevation)
             hs_e = hl
         end if

         if((x.le.0.d0).and.(x.gt.-rl)) then
             hs_w = hl
             hs_e = 0.d0
         end if

         if(x.le.-rl) then
             hs_w = hl
             hs_e = hl - rl*dtan(elevation)
         end if

!- - - - - the hight of ( south, north )'s shade - - -

         if((y.ge.0.d0).and.(y.lt.rl)) then
             hs_s = 0.d0
             hs_n = hl
         end if

         if(y.ge.rl) then
             hs_s = hl - rl*dtan(elevation)
             hs_n = hl
         end if

         if((y.le.0.d0).and.(y.gt.-rl)) then
             hs_s = hl
             hs_n = 0.d0
         end if

         if(y.le.-rl) then
             hs_s = hl
             hs_n = hl - rl*dtan(elevation)
         end if

      end if

!- - - - - - - - - - - - - - - - -

      if(elevation.le.0.d0) then

         ws_e = wl
         ws_w = wl
         ws_s = wl
         ws_n = wl

        else

!- - - - wall east, west - - - - -

         if((x.ge.0.d0).and.(x.lt.rl)) then
             ws_w = 0.d0
             ws_e = wl
         end if

         if(x.ge.rl) then
             ws_w = max(wl-dabs(y),0.d0)
             ws_e = wl
         end if

         if(dabs(y).ge.(wl+rl)) then
             ws_w = wl
         end if

         if((x.le.0.d0).and.(x.gt.-rl)) then
             ws_w = wl
             ws_e = 0.d0
         end if

         if(x.le.-rl) then
             ws_w = wl
             ws_e = max(wl-dabs(y),0.d0)
         end if

         if(dabs(y).ge.(wl+rl)) then
             ws_e = wl
         end if

!- - - - - wall south, north - - - - -

         if((y.ge.0.d0).and.(y.lt.rl)) then
            ws_s = 0.d0
            ws_n = wl
         end if

         if(y.ge.rl) then
            ws_s = max(wl-dabs(x),0.d0)
            ws_n = wl
         end if

         if(dabs(x).ge.(wl+rl)) then
            ws_s = wl
         end if

         if((y.le.0.d0).and.(y.gt.-rl)) then
            ws_s = wl

            ws_n = 0.d0
         end if

         if(y.le.-rl) then
            ws_s = wl
            ws_n = max(wl-dabs(x),0.d0)
         end if

         if(dabs(x).ge.(wl+rl)) then
            ws_n = wl
         end if

     end if

        shade_e = hs_e*wl             !s_e
        shade_w = hs_w*wl             !s_w
        shade_s = hs_s*wl             !s_s
        shade_n = hs_n*wl             !s_n

        sq_wall = hl * wl



        ratio_e = 1.d0 - shade_e/sq_wall  ; ratio(1)=ratio_e
        ratio_w = 1.d0 - shade_w/sq_wall  ; ratio(2)=ratio_w
        ratio_s = 1.d0 - shade_s/sq_wall  ; ratio(3)=ratio_s
        ratio_n = 1.d0 - shade_n/sq_wall  ; ratio(4)=ratio_n


        if(dabs(x).ge.rl) then
           xs = rl
         else

           xs = x
        end if

        if(dabs(y).ge.rl) then
           ys = rl
         else
           ys = y
         end if

     sq_road = 2.d0*rl*wl + rl**2
     sq_shade = wl * ( dabs(xs) + dabs(ys) )

     if(sq_shade>=sq_road) sq_shade = sq_road

     ratio_r  = 1.d0 - sq_shade/sq_road  ; ratio(5) = ratio_r

     endif


!- for output
     if(elevation>=0.d0) then
        write(116,'(9(2x,f13.4))') &
                                dfloat(istep)*dt/3600.d0, x,y,sq_shade, &
                                 ratio(1:5)
     end if

  return
  end subroutine shade
!==============================================================================






















!=========================================================================
  subroutine shade_anal(                                           &
                          istep,                                   & !- in
                          z, dz, hl, wl, rl,                       & !- in
                          elevation, alpha_s, alpha_c,             & !- in
                          ratio_g, ratio                           & !- out
                         )
!=========================================================================
!  to calculate ratio of area under sunlight
!     ratio_g (scalar): for road
!     ratio (direction, canop_layer) : for walls ( 4 direction, kce-2 )
!-------------------------------------------------------------------------

  implicit none

   integer, intent(in)                       :: istep
   real(8), dimension(kme+1), intent(in)     :: z
   real(8), dimension(kme)                   :: dz
   real(8), intent(in)                       :: hl, wl, rl
   real(8), intent(in)                       :: elevation, alpha_s, alpha_c
   real(8), dimension(4,kce-2), intent(out)  :: ratio
   real(8)                                   :: ratio_g
   real(8), dimension(kce-2)                 :: ratio_e, ratio_w, ratio_s, ratio_n

!- local variables

   integer                       :: i, k, id
   real(8)                       :: x, y, xa, ya
   real(8)                       :: xs, ys
   real(8)                       :: sq_shade, sq_road, sq_layer
   real(8)                       :: sq_sun
   real(8)                       :: alpha_tan, alpha_cot, alpha_sin, alpha_cos
   real(8)                       :: a0
   real(8), dimension(4)         :: rate, sq_kabe
   real(8), dimension(kce-2)     :: dzcap, zcap, hs
   real(8), dimension(4,kce-2)   :: sq_wall
   real(8)                       :: tanel
   real(8)                       :: sq_wall_sn1, sq_wall_ew1, h_wall_sn1, h_wall_ew1
   real(8)                       :: sq_wall_sn2, sq_wall_ew2, h_wall_sn2, h_wall_ew2
   real(8)                       :: xx, yy
!--------------------------------------------------------------------------

       do k=1, kce-2
          dzcap(k) = dz(k+1)
          zcap(k)  = z(k+1)
       end do
          sq_kabe(:)   = 0.d0 !hl*wl
          ratio(:,:)   = 0.d0
          rate(:)      = 0.d0
          sq_wall(:,:) = 0.d0
          sq_layer     = wl*dzcap(1)
          sq_road = 2.d0*wl*rl + rl**2
       tanel    = dtan(elevation)  ; if(tanel<=0.1d0) tanel = 0.1d0


!- shadow area
        x =  hl*alpha_s/tanel
        y =  -hl*alpha_c/tanel


        xa = dabs(x)
        ya = dabs(y)

!if( xa >= wl) xa = wl
        if( xa <= 0.001d0) xa = 0.001d0
!if( ya >= wl) ya = wl
        if( ya <= 0.001d0) ya = 0.001d0

          alpha_tan =  xa/ya
          alpha_cot =  ya/xa

          a0 = ( xa + ya ) * wl

        if( x < 0.d0 .and. y < 0.d0)   id = 1  ! shade in
        if( x < 0.d0 .and. y >= 0.d0)  id = 2
        if( x >= 0.d0 .and. y >= 0.d0) id = 3
        if( x >= 0.d0 .and. y < 0.d0)  id = 4



!-----------------------------------------------------
!-               before sun rise
!-----------------------------------------------------
  if ( elevation.lt.0.0d0 ) then
                ratio(:,:) = 0.d0
                ratio_g    = 0.0d0
    else


! when shadow does not reach the other
! building---------------------------
          if ( xa <= rl .and. ya <= rl ) then

              sq_shade   = a0

            else if ( xa > rl .and. ya <= rl) then

              sq_shade    =  a0 - (xa - rl)*(wl-rl*alpha_cot)

            else if ( xa <= rl .and. ya > rl) then

              sq_shade    =  a0 - (ya - rl) * (wl-rl*alpha_tan)

! when the shadow reaches the other three buildings
            else if ( xa > rl .and. ya > rl) then

              sq_shade = 300.d0
!a0 - (xa - rl) * (wl-rl*alpha_cot)  &
!   - (ya - rl) * (wl-rl*alpha_tan)  &
!   - (xa - rl)*(ya -rl)

          end if

             if( sq_shade >= sq_road ) sq_shade = sq_road
             if( sq_shade <= 0.d0) sq_shade = 0.d0
             ratio_g  = 1.d0 - sq_shade/sq_road
!-------------------------------------------------------------------------




       do i=1, 2


! when shadow does not reach the other
! building---------------------------
          if ( xa <= rl .and. ya <= rl ) then

              sq_kabe(i)   = 0.d0

            else if ( xa > rl .and. ya <= rl) then

              xx = hl /xa

              sq_kabe(i)   =  (xa - rl)*(wl-rl*alpha_cot) * xx

            else if ( xa <= rl .and. ya > rl) then
              yy = hl /ya

              sq_kabe(i)   =  0.d0 !(ya - rl) * (wl-rl*alpha_tan) * yy

! when the shadow reaches the other three buildings
            else if ( xa > rl .and. ya > rl) then
              xx = hl /xa
              sq_kabe(i) = ( (xa - rl) * (wl-rl*alpha_cot) * xx  &
                            + (xa - rl)*(ya -rl) ) * xx*0.5d0

          end if

!if( sq_kabe(i) >= wl*hl ) sq_kabe(i) = wl*hl
             if( sq_kabe(i) <= 0.d0) sq_kabe(i) = 0.d0

             rate(i) = 1.d0 - sq_kabe(i)  / (hl*wl)

       end do





       do i=3, 4


! when shadow does not reach the other
! building---------------------------
          if ( xa <= rl .and. ya <= rl ) then

              sq_kabe(i)   = 0.d0

            else if ( xa > rl .and. ya <= rl) then

              xx = hl /xa

              sq_kabe(i)   =  0.d0  !(xa - rl)*(wl-rl*alpha_cot) * xx

            else if ( xa <= rl .and. ya > rl) then
              yy = hl /ya

              sq_kabe(i)   = (ya - rl) * (wl-rl*alpha_tan) * yy

! when the shadow reaches the other three buildings
            else if ( xa > rl .and. ya > rl) then
              yy = hl /ya
              sq_kabe(i) =  ( (ya - rl) * (wl-rl*alpha_tan) * yy &
                            + (xa - rl)*(ya -rl) ) * yy*0.5d0

          end if

!if( sq_kabe(i) >= wl*hl ) sq_kabe(i) = wl*hl
             if( sq_kabe(i) <= 0.d0) sq_kabe(i) = 0.d0

             rate(i) = 1.d0 - sq_kabe(i)  / (hl*wl)

       end do






         if( id == 1 ) then     ! sun in north-east
            rate(2) = 0.d0      ! west
            rate(3) = 0.d0      ! south

         end if
         if( id == 2 ) then     ! sun in south-east
            rate(2) = 0.d0      ! west
            rate(4) = 0.d0      ! north
         end if
         if( id == 3 ) then     ! sun in south-west
            rate(1) = 0.d0      ! east
            rate(4) = 0.d0      ! north
         end if
         if( id == 4 ) then     ! sun in north-west
            rate(1) = 0.d0      ! east
            rate(3) = 0.d0      ! south
         end if







       do k=1, kce-2
         if( id == 1 ) then     ! sun in north-east
           ratio(2,k) = 0.d0
           ratio(3,k) = 0.d0

         end if
         if( id == 2 ) then     ! sun in south-east
           ratio(2,k) = 0.d0
           ratio(1,k) = 0.d0
         end if
         if( id == 3 ) then     ! sun in south-west
           ratio(1,k) = 0.d0
           ratio(4,k) = 0.d0
         end if
         if( id == 4 ) then     ! sun in north-west
           ratio(1,k) = 0.d0
           ratio(2,k) = 0.d0
         end if


       end do



!-------
  end if
!-------



end subroutine shade_anal



!===============================
  end module module_sf_driver
!===============================
