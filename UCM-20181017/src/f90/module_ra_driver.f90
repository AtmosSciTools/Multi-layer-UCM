module module_ra_driver
!=====================================
use module_ra_kondo94
use module_params, only :             &
              kms, kme, lon, lat, dt
use module_io
!===========
contains
!===========
!=============================================================================
subroutine radiation_driver (                                         &
                              istep, z, dz, dz_t,                     & !- in
                              year, month, day, hour, minute, second, & !- in
                              p, rho_a, t, qv,                        & !- in
                              adi, alpha_s, alpha_c, elevation,       & !- inou
                              sw, sd, ss, lw,                         & !- inou
                              sw_ra_heating, lw_ra_heating            & !- inou
                            )
!  to calculate solar elevation, solar azimuth
!               solar short wave radiation ( direct and diffuse on horizon)
!               solar long wave radiation
!               in this model, sw_ra_heating and lw_ra_heating was not
!                   calculate  ( 0 as default)
!.............................................................................
implicit none
!- input variables
integer, intent(in)                             :: istep
integer, intent(in)                :: year, month, day, hour, minute, second
real(8), dimension(kms:kme), intent(in)       :: dz, dz_t
real(8), dimension(kms:kme+1), intent(in)     :: z
real(8), dimension(kms:kme), intent(in)       :: rho_a   ! air density (kg/m3)
real(8), dimension(kms:kme), intent(in)       :: t, p
real(8), dimension(kms:kme), intent(in)       :: qv
!- output variables
real(8), intent(inout)                        :: sw, sd, ss, lw
real(8), intent(inout)                        :: adi, alpha_s, alpha_c, &
                                                   elevation
real(8), dimension(kms:kme), intent(inout)    :: sw_ra_heating, lw_ra_heating
!- local variables
real(8)                            :: humid, albedo, bdust
real(8)                            :: s0
real(8)                            :: xx, lat_local, lon_local, hour_angle
!-------------------------------------------------------------------------------
xx = elevation
albedo = 0.2d0
!------------------------------------
! Calculate radiation flux at top
!------------------------------------
call cal_sw_top (                                          &
                  istep,                                  & !- in
                  year, month, day, hour, minute, second, & !- inout
                  adi, alpha_s, alpha_c, elevation, s0    & !- inout
                 )

bdust  = bdust_urban
humid  = humid_urban
sw = 0.d0
sd = 0.d0
ss = 0.d0
lw = 0.d0

!------------------------------------------
! Calculate radiation flux at the surface
!------------------------------------------
call ra_sw_kondo94 (                        &
                     istep, z, dz,          & !- in
                     s0, elevation,         & !- in
                     humid, albedo, bdust,  & !- int
                     t,                     & !- in
                     sw, sd, ss             & !- inout
                     )

sw_ra_heating(:) = 0.d0

call ra_lw_kondo94 (                           &
                     istep, z, dz, rho_a,      &
                     t, qv, p, humid,          & !- in
                     lw                        & !- inout
                   )

!lw = 300.0d0
lw = lw*0.85
sw = sw*0.85
sd = sd*0.85
ss = ss*0.85
lw_ra_heating(:) = 0.d0

return
end subroutine radiation_driver
!==============================================================================












!============================================================================
  subroutine cal_sw_top (                                         &
                          istep,                                  & !- in
                          year, month, day, hour, minute, second, & !- inout
                          adi, alpha_s, alpha_c, elevation, s0    & !- out
                        )
!============================================================================
!  this program to calc. s0: short wave radiation on horizontal surface on top
!                        elevation: solar elevation angle ( 0 <  < 90 deg)
!                        adi: solar azimuth angle
!                        alpha_s, alpha_c is sin and cosin of solar azimuth
!                        angle
!----------------------------------------------------------------------------

  implicit none

   integer, intent(in)       :: istep
   integer, intent(in)       :: year, month, day, hour, minute, second
   real(8), intent(inout)    :: s0
   real(8), intent(inout)    :: elevation
   real(8), intent(inout)    :: adi
   real(8), intent(inout)    :: alpha_s, alpha_c

!- local variables

   integer :: im
   integer :: mday1(12), mday2(12)
   integer :: leap_year
   real(8) :: doy
   real(8) :: alat
   real(8) :: pi
   real(8) :: theta
   real(8) :: solar_decline
   real(8) :: d1
   real(8) :: et
   real(8) :: hour_angle
   real(8) :: elevation_s
   real(8) :: gmt_angle

   data mday1 / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
   data mday2 / 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

!..........................................................................


       pi = 4.d0 *datan(1.d0)

! leap year
       leap_year = 0
       if(mod(year,4).eq.0) leap_year = 1
       if(mod(year,100).eq.0) leap_year = 0
       if(mod(year,400).eq.0) leap_year = 1

!- return days of month to days of year (sekisan nissu)
       call cal_doy (                                              &
                      leap_year, month, day, hour, minute, second, & !- in
                      doy                                          & !- out
                     )

!- equation of time (kinzisa)
      theta = 2.d0*pi/365.d0 * ( doy - 1.d0 )

      if(leap_year.eq.1) theta = 2.d0*pi/366.d0 * ( doy - 1.d0 )

!- declination of the sun ( do lech)
      solar_decline = 0.006918d0 - 0.399912d0*dcos(theta) + 0.070257d0*dsin(theta) &
                     - 0.006758d0*dcos(2.d0*theta) + 0.000907d0*dsin(2.d0*theta)  &
                     - 0.002697d0*dcos(3.d0*theta) + 0.001480d0*dsin(3.d0*theta)


!-(d0/d1)
      d1 = 1.000110d0 + 0.034221d0*dcos(theta) + 0.001280d0*dsin(theta) &
            + 0.000719d0*dcos(2.d0*theta) + 0.000077d0*dsin(2.d0*theta)

      et = 0.000075d0 + 0.001868d0*dcos(theta) - 0.032077d0*dsin(theta) &
            - 0.014615d0*dcos(2.d0*theta) - 0.040849d0*dsin(2.d0*theta)



!-hour angle (zikaku)
      gmt_angle   = 135.d0
      hour_angle  = ( (dfloat(hour)+dfloat(minute)/60.d0+dfloat(second)/3600.d0-12.d0 ) &
                      + ( lon - gmt_angle )/15.d0 ) * pi/12.d0         !+ et

!- latitude by radian
      alat = pi/180.d0 * lat


! >>>>> ...........................................................<<<<<< !

!- sin of solar elevation angle
      elevation_s =  dsin(alat)*dsin(solar_decline) + dcos(alat)*dcos(solar_decline) * dcos(hour_angle)

!- solar elevation angle
      elevation   =  dasin(elevation_s)

!- solar azimud angle (rad)
      adi = datan(dcos(alat)*dcos(solar_decline)*dsin(hour_angle) &
            / (dsin(alat)*dsin(elevation)-dsin(solar_decline)))


!- solar azimuth angle (sin and cosin)--
!  when alpha_s = plus, solar is in the west, and vice versa
!  when alpha_c = plus, solar is in the north, and vice versa
!  e.g.: in summer morning, solar is in the north-east, then
!        alpha_s is minus, and alpha_c is plus

     alpha_s = dcos(solar_decline)*dsin(hour_angle)/dcos(elevation)

     alpha_c = dcos(alat)*dsin(solar_decline)/dcos(elevation) &
              - dsin(alat)*dcos(solar_decline)*dcos(hour_angle)/dcos(elevation)

     if(alpha_s.ge.1.d0) alpha_s  = 1.d0
     if(alpha_s.le.-1.d0) alpha_s = -1.d0
     if(alpha_c.ge.1.d0) alpha_c  = 1.d0
     if(alpha_c.le.-1.d0) alpha_c = -1.d0


! solar short wave radiation on horizontal surface in the top of air layer.

     s0 = 1365.d0 * d1 * dsin(elevation)

!------------------------------------------------------------
     call put_var_nc(istep, "ELEVATION", elevation)
     call put_var_nc(istep, "ADI", adi)
     call put_var_nc(istep, "COSADI", alpha_c)
     call put_var_nc(istep, "SINADI", alpha_s)



  return
  end subroutine cal_sw_top
!==============================================================================









!=======================================================================
  subroutine cal_doy (                                               &
                        leap_year, month, day, hour, minute, second, & !- in
                        doy                                          & !- out
                     )
!=======================================================================
! subroutine to calculate days of year
!            days of year is accumlated days from 01/01
!            included hour, min, sec information as fraction of day
!            leap year was considered.
!-----------------------------------------------------------------------

 implicit none

   integer, intent(in)    :: month, day, hour, minute, second
   integer, intent(in)    :: leap_year
   real(8), intent(out)   :: doy

!-- local variables

   integer    :: im
   integer    :: mday(12)
   data mday / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
!.......................................................................


   doy = 0.d0

   if(month.gt.1) then
    do im = 1,month-1
     doy = doy + mday(im)
     if(im.eq.2) doy = doy + leap_year
    enddo
   endif

   doy = doy + dfloat(day) + dfloat(hour)/24.d0           &
         + dfloat(minute)/1440.d0 + dfloat(second)/86400.d0


  return
  end subroutine cal_doy

!==============================================================================




end module module_ra_driver
