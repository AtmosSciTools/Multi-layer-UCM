!=============================
module module_pbl_driver
contains
!==================================================
subroutine  pbl_driver(                   &
                       istep,             &  !-in
                       z, dz, dz_t,       &  !-in
                       u, v, t,           &  !-in
                       km, kh             &  !-out
                      )
!==================================================
use module_params,only: kms, kme, bl_pbl_physics 
implicit none
integer,intent(in)                      ::  istep
real(8),dimension(kms:kme),intent(in)   ::  dz, dz_t
real(8),dimension(kms:kme),intent(in)   ::  u, v, t
real(8),dimension(kms:kme+1),intent(in) ::  z
real(8),dimension(kms:kme),intent(out)  ::  km, kh
!---------------------------------------------------
pbl_scheme: select case( bl_pbl_physics )
case(0)
   call pbl_const( km, kh )
case(1) 
   call pbl_my_lev2  (                     &
	               istep,              & !- in
	               z, dz, dz_t,        & !- in
	               u, v, t,            & !- in
	               km, kh              & !- out
                      )
case(2)    
end select pbl_scheme 
return
end subroutine pbl_driver
!====================================================














!=======================================================
subroutine pbl_const( km, kh )  !- inout
use module_params,only:   kms,kme
implicit none
real(8), dimension(kms:kme), intent(inout)   :: kh, km
!- local variables
integer ::  k
do k = kms, kme
   km(k) = 5.0d0
   kh(k) = 5.0d0 
enddo
return
end subroutine pbl_const
!=======================================================







!==========================================================
     subroutine  pbl_test(                       &
                              istep,             &  !-in
                              z, dz, dz_t,       &  !-in
                              u, v, t,           &  !-in
                              km, kh             &  !-out
                             )
!==========================================================

 use module_params,only:                                 &
                        kms,kme, g, kappa, tb,           &
                        ganma, history_interval


   implicit none

    integer,intent(in)                      ::  istep
    real(8),dimension(kms:kme),intent(in)   ::  dz, dz_t
    real(8),dimension(kms:kme),intent(in)   ::  u, v, t
    real(8),dimension(kms:kme+1),intent(in) ::  z
    real(8),dimension(kms:kme),intent(out)  ::  km, kh
 

!- local variables

    integer                            ::   k
    real(8),dimension(kms:kme)         ::   dif1_u, dif1_v,dif1_t, gam
    real(8)                            ::   ref, difmax, difmin

    real(8), parameter  ::   a1 = 0.78d0
    real(8), parameter  ::   a2 = 0.78d0
    real(8), parameter  ::   b1 = 15.d0
    real(8), parameter  ::   b2 = 8.d0
    real(8), parameter  ::   c1 = 0.056d0
    parameter( difmax = 300.d0, difmin = 1.0d0 )
 
    real(8)             ::   x1,   x2,  x3,  x4
    real(8)             ::   gm1, gm2 
    real(8)             ::   smm,  sm, shh, sh
    real(8)             ::   theta0
    real(8)             ::   mls                       !- master length scale

    real(8)             ::   rff, rf, ri
!------------------------------------------------------------------------------


!- and set variables

	   gm1     =      1.d0/3.d0  -  ( 2.d0 * a1/b1 )
	   gm2     =      b2 / b1    +   6.d0  * a1/b1
	   theta0  =      tb


    do k = kms, kme-1

!- master length scale
	     mls =    kappa * z(k+1) /( 1.d0 + kappa * z(k+1) / 100.d0) 
!- Diffusion term calculate
	     dif1_u(k)  =  ( u(k+1) -  u(k) ) / dz_t(k)
	     dif1_v(k)  =  ( v(k+1) -  v(k) ) / dz_t(k)
	     dif1_t(k)  =  ganma  + ( t(k+1) -  t(k) ) / (dz_t(k))

	     ref    =  dif1_u(k)**2d0  +  dif1_v(k)**2d0

!-Data Ajustment-----------------------
       if(ref<=1.d-20) ref   =  1.d-20      
!--------------------------------------


         ri   =  (g/theta0) * dif1_t(k)/ref            !-  Richardson number Ri

         rf    =  0.725d0*                              &
              ( ri+0.186d0 - (ri**2d0- 0.316d0*ri + 0.0346d0)**0.5d0 )


!- Data Ajustment by Rf value
            if(rf > 0.213d0) then
                km(k)     =   difmin
                kh(k)     =   difmin
             else
       
!---------------------------
         
         gam(k) =  rf /( 1.d0 - rf )  
 

         x1    =    gm1  -  c1  - ( 6.d0 * a1  +  3.d0 * a2 )* gam(k) / b1 
         x2    =    gm1  -  gm2 * gam(k)  +  3.d0 * a1 * gam(k) / b1 
         x3    =    gm1  -  gm2 * gam(k) 

         smm   =    3.d0 * a1 * x1 * x3/x2 
         shh   =    3.d0 * a2 * ( gm1 - gm2 * gam(k) )



!- Data Ajustment by Smm, Shh
           if(smm <= 1.d-20) smm = 1.d-20                          
           if(shh <= 1.d-20) shh = 1.d-20
!-----------------------------------

        sm    =  b1**0.5d0 * ( 1.d0 - rf )**0.5d0 * smm**1.5d0
           if(sm <= 1.d-20) sm = 1.d-20 
        sh    =  sm * shh / smm
           if(sh <= 1.d-20) sh = 1.d-20

!- Calculate  Diffution Coefficient  K
 
           km(k)   =   mls**2d0 * ref**0.5d0 * sm
           kh(k)   =   mls**2d0 * ref**0.5d0 * sh


!--------
	     end if
!--------

!  Data Ajustment
!- km by km value
             if(km(k)<= difmin)  km(k)  =  difmin
             if(km(k) >=difmax) km(k)  =  difmax
!- kh by itself value
             if(kh(k)<= difmin)  kh(k)  =  difmin
             if(kh(k) >=difmax/0.76d0 ) kh(k)  =  difmax/0.76d0
!--------------------------------------------
   

!
    if(mod(istep,history_interval)==0) then
 
       write(300,'(f4.1,f7.1,f9.5, f12.5, 3f10.5,f10.7)') &
                                  float(istep)/float(history_interval),&
                                  z(k), u(k+1)-u(k), v(k+1)-v(k), sm,  &
                                  rf, &
                                  ganma  + ( t(k+1) - t(k) ) / dz_t(k), &
                                  ref
    
     end if
!

   end do


  return
end subroutine pbl_test
!================================================================












!===========================================================
subroutine pbl_my_lev2  (                     &
                            istep,            & !- in
                            z, dz, dz_t,      & !- in
                            u, v, t,          & !- in
                            km, kh            & !- out
                          )
!===========================================================
use module_params,only:                                 &
                       kms,kme, g, kappa, tb,           &
                       ganma
implicit none
integer, intent(in)                                :: istep
real(8), dimension(kms:kme+1), intent(in)          :: z
real(8), dimension(kms:kme), intent(in)            :: dz, dz_t
real(8), dimension(kms:kme), intent(in)            :: u, v, t
real(8), dimension(kms:kme), intent(inout)         :: km, kh
!- local variables
integer                              :: k
real(8)                              :: l
real(8)                              :: u1, v1
real(8), dimension(kms:kme)          :: dudz, dvdz, dwdz
real(8), dimension(kms:kme)          :: dtdz, dtdz_r
real(8)                              :: dudz_c, dvdz_c, dtdz_c
real(8)                              :: sq_dudz, sq_u
real(8)                              :: ri, rf
real(8)                              :: gan, gan1, gan2
real(8)                              :: smt, sht, sm, sh
real(8)                              :: theta0
real(8)                              :: difmax, difmin
real(8)                              :: rfc
real(8)                              :: a1, a2, b1, b2, c1
parameter( rfc = 0.213d0 )
parameter( difmax = 300.d0, difmin = 0.5d0 )
parameter( a1=0.78d0, a2=0.78d0, b1=15.d0, b2=8.d0, c1=0.056d0 )

theta0   = tb
!--------- Km ------------
do k = kms, kme-1
        
   dudz(k) = ( u(k+1) - u(k) ) / dz_t(k)
   dvdz(k) = ( v(k+1) - v(k) ) / dz_t(k)
   dtdz(k) = ( t(k+1) - t(k) ) / dz_t(k)
   dtdz_r(k) = ( (ganma*(z(k+2)-0.5d0*dz(k+1)) + t(k+1)) &
                - (ganma*(z(k+1)-0.5d0*dz(k)) + t(k)) ) / dz_t(k)
end do
!-
do k = kms, kme-1
   l = kappa * z(k+1) / ( 1.d0 + kappa * z(k+1) / 100.d0 )
   dudz_c = dudz(k) 
   dvdz_c = dvdz(k) 
   dtdz_c = dtdz(k)

   sq_dudz = dudz_c**2 + dvdz_c**2 
   if(sq_dudz.lt.1.d-10) sq_dudz = 1.d-10

   ri = (g/theta0) * dtdz_r(k) / sq_dudz
   rf = 0.725d0 * ( ri+0.186d0 - ( ri**2d0-0.316d0*ri+0.0346d0 )**0.5d0 )

   if(rf.gt.rfc) then

      km(k)  = difmin
      kh(k)  = difmin

   else

      gan = rf / ( 1.d0 - rf )

      gan1 = 1.d0/3.d0 - 2.d0*a1/b1
      gan2 = b2/b1 + 6.d0*a1/b1

      smt = 3.d0*a1 * ( gan1 - gan2*gan ) & 
            * ( gan1 - c1 - ( 6.d0*a1 + 3.d0*a2 ) * ( gan/b1 ) ) &
            / ( gan1 - gan2*gan + 3.d0*a1 * ( gan/b1 ) )
      sht = 3. * a2 * ( gan1 - gan2*gan )

      if(smt.lt.1.d-10) smt = 1.d-10
      if(sht.lt.1.d-10) sht = 1.d-10

      sm = ( b1 * ( 1. - rf ) * smt**3d0 )**0.5d0   
      if(sm.lt.1.d-10) sm = 1.d-10
      sh = sm * sht / smt                           
      if(sh.lt.1.d-10) sh = 1.d-10
         
      km(k) = (l**2d0) * (dudz_c**2d0 + dvdz_c**2d0)**0.5d0 * sm 
      kh(k)  = (l**2d0) * (dudz_c**2d0 + dvdz_c**2d0)**0.5d0 * sh

   endif

   if(km(k).gt.difmax) km(k) = difmax
   if(km(k).lt.difmin) km(k) = difmin

   if(kh(k).gt.(difmax/0.76d0)) kh(k) = difmax/0.76d0
   if(kh(k).lt.difmin) kh(k) = difmin

enddo
return
end subroutine pbl_my_lev2
!============================






end module module_pbl_driver
