module module_phys_funcs
use module_params, only :  &
         history_interval 
use module_vars, only : istep
contains

!------------------------------------------
function q_rh ( rh, es, t )

   implicit none
   real(8)    :: q_rh
   real(8), intent(in)  :: t, rh, es

   q_rh = 217.d0 * es/(t+273.15d0)*rh/100.d0
   end function q_rh
!-------------------------------------------




!============================================================
SUBROUTINE mos( B1, Z, Z0, UA, TA, TSF,  & ! in
               RIB, CD, CH, XXX )      ! OUT

!  XXX:   z/L (requires iteration by Newton-Rapson method)
!  B1:    Stanton number
!  PSIM:  = PSIX of LSM
!  PSIH:  = PSIT of LSM
IMPLICIT NONE
REAL(8), INTENT(IN)    :: B1, Z, Z0, UA, TA, TSF 
REAL(8), INTENT(OUT)   :: CD, CH
REAL(8), INTENT(INOUT) :: XXX, RIB
REAL(8)                :: XXX0, X, X0, FAIH, DPSIM, DPSIH
REAL(8)       :: F, DF, XXXP, US, TS, AL, XKB, DD, PSIM, PSIH
REAL(8)                :: ALPHA
INTEGER                :: NEWT
INTEGER, PARAMETER     :: NEWT_END=10

IF(RIB <= -15.D0) RIB=-15.D0 
! in the case rib is minus, i.e., unstable
IF(RIB < 0.) THEN
   DO NEWT=1,NEWT_END
      IF(XXX >= 0.) XXX=-1.E-3

      XXX0=XXX*Z0/(Z+Z0)

      X=(1.-16.D0*XXX)**0.25D0
      X0=(1.-16.D0*XXX0)**0.25D0

      PSIM = DLOG((Z+Z0)/Z0)               &
             -DLOG((X+1.)**2.*(X**2.+1.))  &
             +2.*ATAN(X)                   &
             +DLOG((X+1.)**2.*(X0**2.+1.)) &
            -2.*ATAN(X0)
      FAIH=1./SQRT(1.-16.*XXX)
      PSIH=DLOG((Z+Z0)/Z0)+0.4*B1          &
          -2.*DLOG(SQRT(1.-16.*XXX)+1.)    &
          +2.*DLOG(SQRT(1.-16.*XXX0)+1.)

      DPSIM=(1.-16.*XXX)**(-0.25)/XXX     &
            -(1.-16.*XXX0)**(-0.25)/XXX
      DPSIH=1./SQRT(1.-16.*XXX)/XXX       &
            -1./SQRT(1.-16.*XXX0)/XXX

      F=RIB*PSIM**2./PSIH-XXX

      DF=RIB*(2.*DPSIM*PSIM*PSIH-DPSIH*PSIM**2.) &
          /PSIH**2.-1.

      XXXP=XXX
      XXX=XXXP-F/DF
      IF(XXX <= -10.) XXX=-10.

   END DO

ELSE IF(RIB >= 0.142857) THEN

   XXX=0.714
   PSIM=DLOG((Z+Z0)/Z0)+7.*XXX
   PSIH=PSIM+0.4*B1

ELSE   ! in the case 0. < rib < 0.142857

   AL=DLOG((Z+Z0)/Z0)
   XKB=0.4*B1
   DD=-4.*RIB*7.*XKB*AL+(AL+XKB)**2.
   IF(DD <= 0.) DD=0.
   XXX=(AL+XKB-2.*RIB*7.*AL-SQRT(DD))/(2.*(RIB*7.**2-7.))
   PSIM=DLOG((Z+Z0)/Z0)+7.*MIN(XXX,0.714)
   PSIH=PSIM+0.4*B1

END IF


US=0.4*UA/PSIM             ! u*
IF(US <= 0.01) US=0.01
TS=0.4*(TA-TSF)/PSIH       ! T*
CD=US*US/UA**2.            ! CD
!ALPHA=RHO*CP*0.4*US/PSIH      ! RHO*CP*CH*U
!ALPHA = RHO*0.4*US/PSIH       ! RHO*CP*CH*U
!CH   = 0.4*US/PSIH
!CH   =  US*TS / UA / (TA-TSF)
CH   =  0.4d0*US/PSIH / UA 
!CH   =  US/PSIH / UA

return 
END SUBROUTINE mos
!=========================================================














!=========================================================
  subroutine louis79_scheme (                    &
                               rib, z_rib, z0,   & !- in
                               ch, cm            & !- out
                            )

  implicit none

   real(8), intent(in)            :: rib, z_rib, z0
   real(8), intent(out)           :: ch, cm

!- local variables
   integer                        :: i
   real(8)                        :: a2
   real(8)                        :: xx, cmb, chb
   real(8), PARAMETER             :: ch_max = 0.1d0
   real(8), PARAMETER             :: cm_max = 0.1d0
!---------------------------------------------------------
! to calculate exchange coefficients using bulk richardson number (BRN)
! roughness length, and height of reference point.
!

    a2 = (0.4d0/dlog(z_rib/z0))**2
!call check_value ( istep, rib )

! for BRN larger than 0
    if(rib.ge.0.d0) then
! Furthermore, for BRN larger than 0.143
       if(rib.ge.0.142857d0) then

          xx = 0.714d0

       else

          xx = rib*dlog( z_rib / z0 )/(1.d0-7.d0*rib)
      

       endif
    
        ch = 0.16d0/0.74d0/  &
             (dlog( z_rib/z0 )+7.d0*min(xx,0.714d0))**2
        cm = 0.16d0/  &
             (dlog( z_rib/z0 )+7.d0*min(xx,0.714d0))**2

    else

        cmb = 7.4d0 * a2 * 9.4d0 * dsqrt( z_rib/z0 )
        chb = 5.3d0 * a2 * 9.4d0 * dsqrt( z_rib/z0 )

        ch  = a2 / 0.74d0 *(1.d0-9.4d0*rib /(1.+chb*dsqrt(-rib)))
        cm  = a2 * (1.d0-9.4d0*rib/(1.d0+chb*dsqrt(-rib)))

    endif
   
    call check_value ( istep, 'check ch', ch )
    call check_value ( istep, 'check cm', cm )


    if( ch >= ch_max ) ch = ch_max
    if( cm >= cm_max ) cm = cm_max


   return 
   end subroutine louis79_scheme
!=======================================================================






!============================================
subroutine heat_1d (dx, kw, tw)
use module_params, only : dt, nw
!  solve heat conduction eq. by 1 dimension
implicit none
real(8), dimension(nw), intent(in)     :: dx   
real(8), dimension(nw), intent(inout)  :: tw, kw
!- local vars
real(8),dimension(nw) :: dxx, a, b, c, d, gg, ss
integer               :: k

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
end subroutine heat_1d
!===============================================================







!=====================================================================

   subroutine sf_temp_newton_bowen ( rdown, ta, tg, hc, gc, bo, x  )

        use module_params, only : tb

        implicit none   
        real(8), intent(in)      :: rdown, ta, tg, hc, gc, bo
        real(8), intent(out)     :: x
! local variables
        integer                  :: i
        real(8)                  :: fx, ffx, fxn, c
        real(8)                  :: a, d, e
        real(8), PARAMETER       :: err = 1.0d-10
        real(8), PARAMETER       :: eps = 0.9d0 
        real(8), PARAMETER       :: sigma = 5.67d-8
        REAL(8)                  :: RR, HR, ELER, GOR, F
        REAL(8)                  :: DRR, DHR, DELER, DGOR, DF, QS

         a =  - eps * sigma 
         d =  - (1.0d0 + 1/bo ) * hc   -  gc 
         c =    gc * ( tg + tb )   

!- initial guess
       x = 300.0d0

!- loop

       do i=1, 100
! chuan bi
         RR   =  RDOWN - EPS*SIGMA*X**4
         HR   =  HC  * ( X  - ( ta  + tb  ) )
         ELER =  (1.d0/bo) * HR 
         GOR  =  GC  * ( X  - (TG+TB) )
 
         DRR  =  -4.D0 * EPS*SIGMA*X**3
         DHR  =  HC
         DELER = HC * (1.d0/bo) 
         DGOR = GC

         F    =  RR - HR - ELER - GOR          
         DF   =  DRR  - DHR  - DELER - DGOR  
         X    =  X - F / DF
         if ( abs(F) < err) exit

      end do     
          


      return

   end subroutine sf_temp_newton_bowen
!=================================================================



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sf_temp_newton_bowen_ad (                     &
                                rdown, eps, ta, tg,      & ! in
                                hc, gc, bo, tb,          & ! in
                                ts, rnet, h, le, g0      & ! inout
                                   )
! this subroutine is to find the surface temperature by solving 
! newton equation
! f(x) = 0 
! f(x) = rnet - h - le - g = 0
! rnet = rdown - eps * sigma * x^4
! h    = hc * (x - ta)
! le   = 1/bo * hr
! g    = gc * (x - tg)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none   
real(8), intent(in)      :: rdown,  &
                            eps,    &
                            ta, tg, &
                            hc, gc, &
                            bo, tb
real(8), intent(inout)   :: ts, rnet, h, le, g0
! local variables
integer                  :: i
real(8), parameter       :: err   = 1.0d-5
real(8), parameter       :: sigma = 5.67d-8
real(8)                  :: x
real(8)                  :: rr, hr, eler, gor, f
real(8)                  :: drr, dhr, deler, dgor, df
!logical                  :: success
!------------------------------------------------------------------



! first guess
x = 0.d0
! loop for newton solving
do i = 1, 100

   rr    =  rdown  -  eps * sigma * (x+tb) **4
   hr    =  hc  * ( x  -  ta  )
   eler  =  (1.d0/bo) * hr 
   gor   =  gc  * ( x  -  tg  )
 
   drr   =  -4.D0 * eps * sigma * (x+tb)**3
   dhr   =  hc
   deler =  hc * (1.d0/bo) 
   dgor  =  gc

   f     =  rr - hr - eler - gor          
   df    =  drr  - dhr  - deler - dgor  
   
   x    =  x - f / df
   
   if ( abs(f) < err) then
       ts = x
       exit
       
   end if

end do
! end loop     
          
rnet    =  rdown  -  eps * sigma * (ts+tb) **4
h       =  hc  * ( ts  -  ta  )
le      =  (1.d0/bo) * hr 
g0      =  gc  * ( ts  - tg )


return

end subroutine sf_temp_newton_bowen_ad
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++














!=====================================================================
subroutine sf_temp_newton_bowen1 (                          &
                                    RDOWN, EPS, TA, TG,    & ! in
                                    HC, GC, LEC, QVA, P,   & ! in
                                    TS, RNET, H, LE, G0    & ! out
                                   )

        use module_params, only : tb
!-----------------------------------------
    implicit none   
       real(8), intent(in)      :: RDOWN,  &  ! downward short wave (W/m2)
                                   EPS,    &  ! emiss
                                   TA,     &  ! air temperature (-tb)
                                   TG,     &  ! Under ground temp. (-tb)
                                   HC,     &  ! (rho x cp x ua x ch)
                                   GC,     &  ! Ground heat exchange co
                                   LEC,    &  ! (rho x le x ua x ch x beta)
                                   QVA,    &  ! air specific humidit (kg/kh)
                                   P         ! pressure (hpa)

       real(8), intent(inout)   :: TS,     &  ! surface temp. (-tb)
                                   RNET,   &  ! net radiation
                                   H,      &  ! sensible heat
                                   LE,     &  ! latent heat
                                   G0        ! ground heat

! local variables
       integer                  :: i
       real(8), PARAMETER       :: ERR   = 1.0d-10
       real(8), PARAMETER       :: SIGMA = 5.67d-8
       REAL(8)                  :: RR, HR, ELER, GOR, F
       REAL(8)                  :: DRR, DHR, DELER, DGOR, DF
       REAL(8)                  :: ES, DES, QS, DQS, X

      
!- initial guess
       x = 300.d0

!- loop
       do i=1, 100
            
             ES  =6.11D0*DEXP( 5404.D0*(X-273.15D0)/(273.15D0*X) )
             DES =5404.D0*ES/(X**2.)
             QS  =0.622D0 * ES/P/ ( 1.d0 - 0.378d0 * (ES / P))
! DQS = DES*0.622*P/((P-0.378*ES)**2.)
! QS  = 0.622*ES/(P-0.378*ES)
! DQS = DES*0.622*P/((P-0.378*ES)**2.)
             DQS = 0.622d0 * 5404.d0 * ES/(P*X**2) 

             RR    =  RDOWN - EPS*SIGMA*X**4 
             HR    =  HC  * ( X  - ( TA  + TB ) )
             ELER  =  LEC * ( QS - QVA ) 
             GOR   =  GC  * ( X  - ( TG  + TB ) )
             DRR   =  -4.D0 * EPS*SIGMA*X**3
             DHR   =  HC
             DELER =  LEC * DQS 
             DGOR  =  GC

           F     =  RR  - HR   - ELER  - GOR        
           DF    =  DRR - DHR  - DELER - DGOR 
           X     =  X - F / DF  
           if ( abs(F) < err .AND. ABS(F/DF) < ERR) exit
           if(X < 250.d0) X = 250.d0    
      end do
           
           TS   = X - TB
           H    = HR
           LE   = ELER
           G0   = GOR
           RNET = RR
           call check_value ( istep, 'check ts newton', ts )


      return

   end subroutine sf_temp_newton_bowen1
!=================================================================





!=====================================================================
subroutine sf_temp_newton (                        &
                            RDOWN, EPS, TA, TG,    & ! in
                            HC, GC, LEC, QVA, P,   & ! in
                            TS, RNET, H, LE, G0    & ! out
                          )

use module_params, only : tb
implicit none   
real(8), intent(in)      :: RDOWN,  &  ! downward short wave (W/m2)
                                   EPS,    &  ! emiss
                                   TA,     &  ! air temperature (-tb)
                                   TG,     &  ! Under ground temp. (-tb)
                                   HC,     &  ! (rho x cp x ua x ch)
                                   GC,     &  ! Ground heat exchange co
                                   LEC,    &  ! (rho x le x ua x ch x beta)
                                   QVA,    &  ! air specific humidit (kg/kh)
                                   P         ! pressure (hpa)

real(8), intent(inout)   :: TS,     &  ! surface temp. (-tb)
                                   RNET,   &  ! net radiation
                                   H,      &  ! sensible heat
                                   LE,     &  ! latent heat
                                   G0        ! ground heat
! local variables
integer                  :: i
real(8), PARAMETER       :: ERR   = 1.0d-10
real(8), PARAMETER       :: SIGMA = 5.67d-8
REAL(8)                  :: RR, HR, ELER, GOR, F
REAL(8)                  :: DRR, DHR, DELER, DGOR, DF
REAL(8)                  :: ES, DES, QS, DQS, X

      
!- initial guess
x = 300.d0
!ES  =6.11D0*DEXP( 5404.D0*(X-273.15D0)/(273.15D0*X) )
!DES =5404.D0*ES/(X**2.)
!QS  =0.622D0 * ES/P/ ( 1.d0 - 0.378d0 * (ES / p))
!DQS = DES*0.622*P/((P-0.378*ES)**2.)
!ES  = 6.11*EXP( (2.5*10.**6./461.51)*(X-273.15)/(273.15*X) )
!DES = (2.5*10.**6./461.51)*ES/(X**2.)
!QS  = 0.622*ES/(P-0.378*ES)
!DQS = DES*0.622*P/((P-0.378*ES)**2.)
!DQS = 0.622d0 * 5404.d0 * ES/(P*X**2)
!call check_value ( istep, 'check es', es )

!- loop
       do i=1, 100
            
             ES  =6.11D0*DEXP( 5404.D0*(X-273.15D0)/(273.15D0*X) )
             DES =5404.D0*ES/(X**2.)
             QS  =0.622D0 * ES/P/ ( 1.d0 - 0.378d0 * (ES / P))
! DQS = DES*0.622*P/((P-0.378*ES)**2.)
! QS  = 0.622*ES/(P-0.378*ES)
! DQS = DES*0.622*P/((P-0.378*ES)**2.)
             DQS = 0.622d0 * 5404.d0 * ES/(P*X**2) 

             RR    =  RDOWN - EPS*SIGMA*X**4 
             HR    =  HC  * ( X  - ( TA  + TB ) )
             ELER  =  LEC * ( QS - QVA ) 
             GOR   =  GC  * ( X  - ( TG  + TB ) )
             DRR   =  -4.D0 * EPS*SIGMA*X**3
             DHR   =  HC
             DELER =  LEC * DQS 
             DGOR  =  GC

           F     =  RR  - HR   - ELER  - GOR        
           DF    =  DRR - DHR  - DELER - DGOR 
           X     =  X - F / DF  
           if ( abs(F) < err .AND. ABS(F/DF) < ERR) exit
           if(X < 250.d0) X = 250.d0    
      end do
           
           TS   = X - TB
           H    = HR
           LE   = ELER
           G0   = GOR
           RNET = RR
           call check_value ( istep, 'check ts newton', ts )

!write(*,'(4f10.3)') x, gor, hr, eler
!if(mod(istep,history_interval*10)==0) then
!  write(*,'(7f15.4)') (X  - ( TA  + TB ))/((QS - QVA)*1000.),  &
!         le, lec, QS - QVA, (X  - ( TA  + TB ))
!end if

!DO ITERATION=1,20
!  ES=6.11*EXP( (2.5*10.**6./461.51)*(TRP-273.15)/(273.15*TRP) )
!  DESDT=(2.5*10.**6./461.51)*ES/(TRP**2.)
!  QS0R=0.622*ES/(PS-0.378*ES)
!  DQS0RDTR = DESDT*0.622*PS/((PS-0.378*ES)**2.)
!  RR=EPSR*(RX-SIG*(TRP**4.)/60.)
!  HR=RHO*CP*CHR*UA*(TRP-TA)*100.
!  ELER=RHO*EL*CHR*UA*BETR*(QS0R-QA)*100.
!  G0R=AKSR*(TRP-TRL(1))/(DZR(1)/2.)
!  F = SR + RR - HR - ELER - G0R
!  DRRDTR = (-4.*EPSR*SIG*TRP**3.)/60.
!  DHRDTR = RHO*CP*CHR*UA*100.
!  DELERDTR = RHO*EL*CHR*UA*BETR*DQS0RDTR*100.
!  DG0RDTR =  2.*AKSR/DZR(1)
!  DFDT = DRRDTR - DHRDTR - DELERDTR - DG0RDTR
!  DTR = F/DFDT
!  TR = TRP - DTR
!  TRP = TR
! IF( ABS(F) < 0.000001 .AND. ABS(DTR) < 0.000001 ) EXIT
! END DO






      return

   end subroutine sf_temp_newton
!=================================================================







!-----------------------------------------------
function qs_dif ( t, p ) 
implicit none
real(8)          :: qs_dif
real(8), intent(in)     :: t, p  
real(8)                 :: a, b, es, es_dif 

! from -30 ~ 100 C
a = 0.622d0 / p
b = - 0.378d0/p
es  = 6.11d0*10.d0**(7.5d0*t/(t+237.3d0)) 

es_dif  =  6.1078d0 * (2500.d0 - 2.4d0 * t) / &
           ( 0.4615d0 * (273.15d0 + t)**2 )   &
           * 10.d0 ** ( 7.5*t / (237.3d0 + t) )

qs_dif = a * es_dif / (1.d0 + b*es)**2
 
end function qs_dif
!------------------------------------------------




!----------------------------------------
function qs_tetens ( t, p )
implicit none
real(8)      :: qs_tetens
real(8), intent(in)   :: t, p
real(8)      :: e_sat
! t: temperature (celsius degree)
! p: air presssure
! e_sat : saturated vapor pressure
e_sat  = 6.11d0*10.d0**(7.5d0*t/(t+237.3d0))
qs_tetens = 0.622d0 * ( e_sat / p) / &
           ( 1.d0 - 0.378d0 * (e_sat / p)) 
end function qs_tetens
!---------------------------------------





















!===================================================================
subroutine beta_method (                      &
                          k_ref,              & ! -in
                          ts, z_t, pres,      & !- in
                          rho, ch, large_u,   & !- in
                          rain, rain_max, qv, & !- in
                          raing,              & !- inout
                          qv_sat, beta,       & !- out
                          le                  & !- out
                         )
use module_params, only : kms, kme, tb, dt
implicit none
integer, intent(in)      :: k_ref
real(8), intent(in)      :: ts, pres, rho, ch, large_u, rain, rain_max
real(8), dimension(kms:kme)   :: qv, z_t
real(8), intent(inout)        :: raing, qv_sat, beta, le
! local variables
real(8)   ::  e_sat, temp, Lv
PARAMETER( Lv = 2.45d6 )         ! Latent heat   [J/kg]
real(8), PARAMETER  :: beta_min = 0.0d0
!-----------------------------------------------------

temp   = ts + tb-273.15d0 - 0.01d0 * z_t(k_ref)
e_sat  = 6.11d0*10.d0**(7.5d0*(temp)/((temp)+237.3d0))
qv_sat = 0.622d0 * ( e_sat / pres) / ( 1.d0 - 0.378d0 * (e_sat / pres)) 

beta   = beta_min

if ( rain > 1.0d0 )  beta = 1.0d0
raing = raing + rain
        
if ( raing >= rain_max ) then
   raing  = rain_max
   beta         = 1.0d0
end if

if ( raing < rain_max ) then
    beta  = raing / rain_max
end if
if ( beta <= beta_min) beta = beta_min 

beta = 0.4d0
le     = Lv * rho * beta * ch * large_u * ( qv_sat - qv(kms+1) )    
! [J/kg x kg/m3 x kg/kg x m/s] = [W/m2]

! kg/m3 x kg/kg x m/s x s = kg/m2 (water/m2)
raing = raing -  rho * beta * ch * large_u * ( qv_sat - qv(kms+1) )  * dt  & 
        * 1.0d0 * 1000.d0             ! kg/m3 x 1000mm  = mm
if(raing<0) raing=0.d0

return
end subroutine beta_method
!==============================================================





!===========================================
subroutine check_value ( istep, mes, x )
implicit none
character (len=*), intent(in)     :: mes
integer, intent(in)     :: istep
real(8), intent(in)     :: x
if(x.ne.x) then  
   write(*,*) istep,  trim(mes), ': it is NaN '      
   stop                                
end if                                 
if((x*0.0d0).ne.0.0d0) then 
   write(*,*) istep, trim(mes), 'it is Inf '      
   stop                                
end if                                 
                
return
end subroutine check_value
!==============================================

end module module_phys_funcs



