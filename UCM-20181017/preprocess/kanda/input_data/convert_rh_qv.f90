program main

implicit none


real(8)     :: rho_a
real(8)     :: p, ta, rh, g
real(8)     :: qv
real(8)     :: z = 250.d0
real(8)     :: sw, sd, ss, lw, ah, u, v , hum



   integer        :: i
   real(8)        :: pres       !- pressure [hPa]
   real(8)        :: e_sat      !- saturated water vapor pressure [hPa]
   real(8)        :: qv_sat     !- saturated mixing ratio of water vapor [kg/kg]
   real(8)        :: temp, theta
   real(8)        :: xx1, xx2
  

  !-- local parameters

   real(8), parameter      :: aa = 7.5d0
   real(8), parameter      :: bb = 237.3d0       !-- for water sureface
   real(8) :: sumqv   

g = 9.8

open(11, file="kanda07.txt")
open(12, file="kanda_data_072930.txt")


do i=1, 97

    read(11,*) sw, sd, ss, lw, ta, hum, u, v, ah,  rh, p 

    rho_a = 1.2d0  
    pres =  p + ( - rho_a * g * z ) * 10.d0**(-2)      ![hPa]
    temp   = ta

                                                                 
    e_sat  = 6.1078d0 * 10.d0 ** (aa* temp /( bb + temp))
   
    qv_sat = 0.622d0 * ( e_sat / pres) / ( 1.d0 - 0.378d0 * (e_sat / pres)) !* 1000.d0
                                                                            ![kg/kg]

    qv  = rh/100.d0 * qv_sat          

  !write(*,'(4f7.3)')  qv,  pres, e_sat,  qv_sat*1000 -hum

   write(12,'(9f10.5)') sw, sd, ss, lw, ta, qv*1000. , u, v, ah  

end do










end program main





