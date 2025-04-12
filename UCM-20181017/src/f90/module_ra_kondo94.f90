module module_ra_kondo94

!=====================================================================
  use module_params, only :                                          &  
    kms, kme, dt, ganma, tb, pb, sigma, g,                           & 
    humid_urban, bdust_urban
!======================================================================

!===========
 contains
!===========






!=================================================================
  subroutine ra_lw_kondo94 (                           &
                             istep, z, dz, rho_a,      &
                             t_1D, qv_1D, p_1D, humid, & !- in
                             lw_0                      & !- inout
                           )
!==================================================================


!------------------------------------------------------------------
!  DISCIPTION
!   Purpose: to calculate
!    from:
!    outside vars: istep: time steps
!                  z, dz: vertical height
!                  rho_a: air humidity
!    outside subs: cal_tdew, cal_wtop
!-----------------------------------------------------------------


  implicit none

   integer, intent(in)                        :: istep
   real(8), dimension(kms:kme), intent(in)    :: dz, rho_a
   real(8), dimension(kms:kme+1), intent(in)  :: z
   real(8), dimension(kms:kme), intent(in)    :: t_1D, qv_1D, p_1D
   real(8), intent(in)                        :: humid
   real(8), intent(inout)                     :: lw_0

!- local variables

   integer  :: i, wt_opt
   real(8)  :: tdew
   real(8)  :: logwt
   real(8)  :: wtop
!----------------------------------------------------------------


  wt_opt = 1

  if ( wt_opt == 1 ) then

    call cal_tdew (                   &
                     dz, t_1D, humid, & !- in
                     tdew             & !- out
                   )

    logwt = 0.0315d0*tdew - 0.1836d0

   else if ( wt_opt == 2 ) then 

     call cal_wtop (                     &
                      z, dz, rho_a,      &
                      t_1D, p_1D, qv_1D, & !- in
                      wtop               & !- out
                   ) 


     logwt = log10(wtop)  
     
   endif

    lw_0 = (0.74d0+0.19d0*logwt+0.07d0*logwt**2)*sigma*tb**4



  return
  end subroutine ra_lw_kondo94
!==============================================================================














!=======================================================
  subroutine ra_sw_kondo94 (                       &
                             istep, z, dz,         & !- in
                             s0, elevation,        & !- in
                             humid, albedo, bdust, & !- int
                             t_1D,                 & !- in
                             sw, sd, ss            & !- inout
                           )
!=======================================================

!-----------------------------------------------------------------!
!  DISCIPTION
!   Purpose: to cal.: sw, sd, ss
!                     sw is global short wave radiation (on horizontal);
!                     sd is direct short wave radiation (vertical);
!                     ss is diffuse short wave radiation (vertical).
!   the amount of solar radiation come to surface, depends on dust
!        vapor, clouds, pollutants in air.
!
!-----------------------------------------------------------------!



  implicit none

!- input variables

   integer, intent(in)                             :: istep
   real(8), intent(in)                             :: s0, elevation      
   real(8), intent(in)                             :: albedo, humid, bdust
   real(8), dimension(kms:kme), intent(in)         :: dz
   real(8), dimension(kms:kme+1), intent(in)       :: z
   real(8), dimension(kms:kme), intent(in)         :: t_1D

!- output variables

   real(8), intent(inout)  :: sw, sd, ss

!- local variables

   real(8)                      :: tdew
   real(8)                      :: p
   real(8)                      :: m, c1, logw
   real(8)                      :: f1, i1, j1
   real(8)                      :: c2, f2, i2

   real(8)                      :: time, sdtop
!-------------------------------------------------------------------------


	    p = pb

	    call cal_tdew (                    &
	                     dz, t_1D, humid,  & !- in
	                     tdew              & !- out
	                   )

! if sun rise
	    if(elevation.gt.0.d0) then

	         m = (p/1000.d0)*(1.d0/dsin(elevation))

	         if(bdust.lt.0.3d0) then
	             c1 = 0.21d0 - 0.2d0 * bdust
	         else
	         c1 = 0.15d0
	         endif

		     logw = 0.0312d0*tdew - 0.963d0
		     f1 = 0.056d0 + 0.16d0 * dsqrt(bdust)
		     i1 = 0.014d0 * (m+7.d0+2.d0*logw) * logw
		     j1 = (0.066d0+0.34d0*dsqrt(bdust))*(albedo-0.15d0)

!- Global solar radiation, after cross air
	        sw = s0 * (c1+0.7d0*10.d0**(-m*f1))*(1.d0-i1)*(1.d0+j1)
	     


	    if(sw.lt.0.d0) sw = 0.d0

	     if(bdust.le.0.3d0) then
	      c2 = 0.15d0 - 0.2d0*bdust
	     else
	      c2 = 0.09d0
	     endif

	     f2 = 0.075d0 + 0.65d0*bdust
	     i2 = 0.02d0*(m+5.5d0+1.5d0*logw)*logw

!- Direct solar radiation
	     sd = s0*(c2+0.75d0*10.d0**(-m*f2))*(1.d0-i2)

	    else

	     sw = 0.d0
	     sd = 0.d0

	    endif
 
             ss = sw - sd
!if(ss<=0.d0 ) ss = 0.d0

  return
 end subroutine ra_sw_kondo94
!==============================================================================







!=======================================================
  subroutine cal_tdew (                   &
                         dz, t_1D, humid, & !- in
                         tdew             & !- out
                       )
!=======================================================

!------------------------------------------------
!   DESCIPTION
!    Purpose: to calculate dew temperature
!     from
!     outside vars; dz: grid distance
!                   t_1D: temperature
!                   humid
!------------------------------------------------


  implicit none

!-  input variables

   real(8), dimension(kms:kme), intent(in) :: dz
   real(8), dimension(kms:kme), intent(in) :: t_1D
   real(8), intent(in)                     :: humid

!- output variables

   real(8), intent(out) :: tdew

!- local variables

   integer                     :: k
   real(8)                     :: gan_t, tc
   real(8)                     :: p
   real(8)                     :: es, e
!--------------------------------------------------------


	   p = pb

	   tc = tb - 273.15d0         !- C degree from basic temperature

	   gan_t = -0.01d0 + ganma    !- ganma of C temp?

	   es = 6.1078d0                                                 &
	       * 10.d0**(7.5d0*(0.5d0*dz(kms+1)*gan_t+t_1D(kms+1)+tc)    &
	            / (237.3d0+(0.5d0*dz(kms+1)*gan_t+t_1D(kms+1)+tc)))  
	   e  = ( humid/100.d0 ) * es
!- dew point temperature
	   tdew = ( 237.3d0*log10(e/6.11d0))/(7.5d0-log10(e/6.11d0))


 return
 end subroutine cal_tdew
!==============================================================================







!=======================================================
 subroutine cal_wtop (                      &
                         z, dz, rho_a,      &
                         t_1D, p_1D, qv_1D, & !- in
                         wtop               & !- out
                       )
!=======================================================

    implicit none

!-  input variables

   real(8), dimension(kms:kme), intent(in) :: dz, rho_a
   real(8), dimension(kms:kme+1), intent(in) :: z
   real(8), dimension(kms:kme), intent(in) :: t_1D, qv_1D, p_1D

!- output variables

   real(8), intent(out) :: wtop  ![cm]

!- local variables

   integer                     :: k
   real(8)                     :: gan_t
   real(8)                     :: p1, p2 
   real(8)                     :: rho_w, wtop2
   real(8)                     :: dr, dr_m, dr_cm
   real(8), dimension(kms:kme) :: press

!------------------------
! Set parameters

	 rho_w = 1000.d0 ![kg/m^3]
	 gan_t = ganma - 0.01d0

!------------------------

	 wtop = 0.d0 
	 wtop2 = 0.d0 
	   
	 do k = kms, kme
	  press(k) = pb*100.d0 + p_1D(k)*rho_a(k) - rho_a(k)*(z(k)+0.5d0*dz(k))
	 end do

	 do k = kms+1, kme-1
	   
	p1 = ( dz(k-1)*press(k) + dz(k)*press(k-1) ) / ( dz(k) + dz(k-1) ) ![Pa]
        p2 = ( dz(k+1)*press(k) + dz(k)*press(k+1) ) / ( dz(k) + dz(k+1) )  ![Pa]
	   dr = - qv_1D(k) * ( p2 - p1 ) / g  !- [kg/m^2]
	   dr_m = dr / rho_w        ![m]
	   dr_cm = dr_m*100.d0      ![cm]
	   wtop = wtop + dr_cm      ![cm]
!   wtop2 = wtop2 + dr*10.d0 ![g/cm^2]

 	end do 


 return
 end subroutine cal_wtop

end module module_ra_kondo94
