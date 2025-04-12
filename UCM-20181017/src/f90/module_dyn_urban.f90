module module_dyn_urban
use module_params, only: ganma, tb, g,  kms, kme,          &
                         ws_opt, t_opt, q_opt, u_opt, v_opt
use module_vars, only : kce,  &
                        THETA_T, QV_T, U_T, V_T, SW_B, LW_B
 contains
!===============================================================
subroutine dyn_urb(                                  &
                    istep,                           & !-in
                    z, z_t, dz, dz_t,                & !-in
                    canop_a,poro, q, qev,            & !-in
                    km, kh,                          & !-in
                    tau_u, tau_v, tau_t, tau_qv,     & !-in
                    u, v, t, qv                      & !-inout
                   )
!===============================================================
! solve dynamic equation of velocity and heat eq.
! use tau, included momentum and heat tau, exhausted from surface
! to calculate u, v, t
! scheme: implicit method for diffusion term
!         explicit method for anothers.
!...............................................................
implicit none
integer,intent(in)                      ::  istep
real(8),dimension(kms:kme),intent(in)   ::  dz, dz_t, z_t
real(8),dimension(kms:kme+1),intent(in) ::  z
real(8),dimension(kms:kme),intent(in)::  canop_a, poro, q, qev
real(8),dimension(kms:kme),intent(in)   ::  km, kh
real(8), intent(in)              ::  tau_u, tau_v, tau_t, tau_qv
real(8),dimension(kms:kme),intent(inout)::  u, v, t, qv
!- local variables
real(8),dimension(kms:kme)        ::  fu, fv, ft, fqv
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!- step one: calculate addition terms
!- in momentum eq: resistance term of building
!- in heat eq: sensible heat flux from building
call bdy_all( istep, z_t, u, v, t, qv)
!call bdy_wrf( istep, z_t, u, v, t, qv)

call dyn_pre (                              &
                dz, dz_t,                   & !- in
                u, v,                       & !- in
                canop_a, q, qev,            & !- in
                fu, fv, ft, fqv             & !- out
              )

!- step two: solve momentum equation
call dyn_momentum(                          &
                dz, dz_t, poro,             & !-in
                km, fu, fv, tau_u, tau_v,   & !-in
                u, v                        & !-inout
                  )

!- solve temperature
call dyn_heat (                             &
                istep, dz, dz_t,            & !- in
                poro, kh, ft, tau_t,        & !- in
                t                           & !- inout
               )

call dyn_moist (                            &
                istep, dz, dz_t,            & !- in
                poro, kh, fqv, tau_qv,      & !- in
                qv                          & !- inout
                )



call bdy_all( istep, z_t, u, v, t, qv )
!call bdy_wrf( istep, z_t, u, v, t, qv)

return
end subroutine dyn_urb
!============================================================





!=====================================================
subroutine bdy_wrf( istep, z_t, u, v, t, qv )

implicit none
integer, intent(in)   :: istep
real(8),dimension(kms:kme), intent(in)   :: z_t
real(8),dimension(kms:kme), intent(inout) :: u, v, t, qv

u(kme)  =  U_T
v(kme)  =  V_T
t(kme)  =  THETA_T + 273.15 - tb   !* GANMA
qv(kme) =  QV_T

return
end subroutine bdy_wrf
!====================================================







!======================================================
subroutine bdy_all( istep,z_t, u, v, t, qv )

use module_params, only:  Ug, Vg
implicit none
integer, intent(in)   :: istep
real(8),dimension(kms:kme), intent(in)    :: z_t
real(8),dimension(kms:kme), intent(inout) :: u, v, t, qv

real(8)      :: tpo
real(8)      :: rho_a
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
! top wind
if ( ws_opt > 1 ) then
   u(kme) = u_t
   v(kme) = v_t
else
   u(kme) = Ug
   v(kme) = Vg
end if

! top theta
rho_a = 1.2d0
if ( t_opt > 1 ) then
   tpo  =  (theta_t + 273.15d0) *   &
           (1000.d0/ (1000.d0   +   &
           ( -rho_a * g * z_t(kme))/100.d0))**(0.286d0)
   t(kme)  =  tpo - tb
else
   t(kme)  =  0.d0  !t(kme-1)
end if

! top hum
if ( q_opt > 1 ) then
   qv(kme) = qv_t
else
   qv(kme) = 0.002d0 !qv(kme-1)
end if
return
end subroutine bdy_all
!======================================================









!===========================================================
 subroutine dyn_pre (                               &
                        dz, dz_t,                   & !- in
                        u, v,                       & !- in
                        canop_a, q, qev,            & !- in
                        fu, fv, ft, fqv             & !- out
                     )
!===========================================================
 use module_params, only:  kms, kme,  dt, cp, lat, omega, &
                           Ug, Vg, sf_surface_physics
!------------------------------------------------------------
!  in momentum equation: > corioris force,
!                          resistance of building ( urban)
!  in heat equation: > heat from building ( urban only )
!------------------------------------------------------------
implicit none
!- inout variables
real(8),dimension(kms:kme),intent(in)    :: dz, dz_t
real(8),dimension(kms:kme),intent(in)    :: u, v
real(8),dimension(kms:kme),intent(in)    :: canop_a, q, qev
real(8),dimension(kms:kme),intent(out)   :: fu, fv, ft, fqv
!- local variables
real(8),dimension(kms:kme)    :: large_u
real(8)                                    :: lat_ra, f
integer                                    :: k
real(8), parameter                         :: cd = 1.d0
!-----------------------------------------------------------

lat_ra  =    3.14159d0 * lat / 180.d0        !  latitude [rad]
f       =    2.d0 * omega * dsin(lat_ra)     !  Corioris [1/s]

if(sf_surface_physics==1) then    ! in case of slab model

  fu(:)  =   f*(v(:)-Vg)      ! only corioris force
  fv(:)  =   -f*(u(:)-Ug)     !   in the momentum eq.
  ft(:)  =   0.0d0            ! no heat frombuilding
  fqv(:) =   0.0d0

else
  large_u(:) = dsqrt(u(:)**2 + v(:)**2)
  where(large_u.lt.1.d0) large_u = 1.d0
! in urban case: fu, fv combine corioris and resistance forces
  fu(:) =   f*(v(:)-Vg)  - cd * canop_a(:) * u(:) *  large_u(:)
  fv(:) =  -f*(u(:)-Ug)  - cd * canop_a(:) * v(:) *  large_u(:)
  ft(:) =   q(:)
  fqv(:) =   0.d0 !qev(k)

end if
return
end subroutine dyn_pre
!==============================================================




!==============================================================
subroutine dyn_heat(                          &
                     istep, dz, dz_t,         & !- in
                     poro, km, ft, tau_t,     & !- in
                     t                        & !- inout
                   )
!==============================================================
! solve heat equation
! use thomat method to solve matrix
use module_params, only:  kms,kme,  dt, cp, sf_surface_physics

implicit none
!- inout variables
integer, intent(in)                        :: istep
real(8),dimension(kms:kme),intent(in)      :: dz, dz_t, km, poro
real(8),dimension(kms:kme),intent(inout)   :: ft
real(8),intent(in)                         :: tau_t
real(8),dimension(kms:kme),intent(inout)   :: t

!- local variables
real(8),dimension(kms:kme) :: a, b, c, d, ntonc
real(8)                    :: x1, x2, x3, x4, tbottom
integer                    :: k
real(8), dimension(:), allocatable    :: aa, bb, cc, dd, xx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!- calculate newton cooling term
ntonc(:) = 0.0d0
ntonc(:) = - t(:)/86400.0d0
!- STEP 1: calculate a, b, c, d coefficients
do k = kms+1, kme-1

   x1 = 0.5d0 * (poro(k-1)+poro(k)) /poro(k)
   x2 = 0.5d0 * (poro(k+1)+poro(k)) /poro(k)
   x4 = (poro(k) - poro(k-1))/ poro(k)
   x3 =  ganma / dz(k) * ( ( km(k) - km(k-1) ) + km(k-1) * x4)

   a(k)   =  ( - km(k-1) / ( dz_t(k-1) * dz(k) ) ) * dt *x1
   c(k)   =  ( - km(k)   / ( dz_t(k)   * dz(k) ) ) * dt *x2
   b(k)   =  1.d0 - a(k)  - c(k)
   d(k)   =  t(k)   + dt * ( ft(k) + x3 + ntonc(k) )

end do

!- boundary condition
!- b.c at bottom

!case of tau from bottom
d(kms+1)   =  d(kms+1)  -   tau_t * (dt/ dz(kms+1))
b(kms+1)   =  1.d0 - c(kms+1)
a(kms+1)   =  0.d0

!- no flux from kme
d(kme-1)     =  d(kme-1) - c(kme-1) * t(kme)
c(kme-1)     =  0.d0
!- STEP 2: solve simultaneous equations by Thormat method
!- solve tridiag
! call solv_tri( a, b, c, d, t)

allocate( aa(kms+1:kme-1))
allocate( bb(kms+1:kme-1))
allocate( cc(kms+1:kme-1))
allocate( dd(kms+1:kme-1))
allocate( xx(kms+1:kme-1))

do k=kms+1, kme-1
   aa(k) = a(k) ; bb(k)=b(k) ; cc(k)=c(k); dd(k)=d(k)
end do

call solv_thomas ( kms+1, kme-1, aa, bb, cc, dd, xx)


t(kms+1:kme-1) = xx(kms+1:kme-1)



deallocate( aa )
deallocate( bb )
deallocate( cc )
deallocate( dd )
deallocate( xx )


return
end subroutine dyn_heat
!================================================================








!================================================================
subroutine dyn_moist(                              &
                     istep, dz, dz_t,              & !- in
                     poro, km, ft, tau_t,          & !- in
                     t                             & !- inout
                     )
!================================================================
! solve heat equation
! use thomat method to solve matrix
use module_params, only:  dt, cp, sf_surface_physics
implicit none
!- inout variables
integer, intent(in)                        :: istep
real(8),dimension(kms:kme),intent(in)      :: dz, dz_t, km, poro
real(8),dimension(kms:kme),intent(inout)   :: ft
real(8),intent(in)                         :: tau_t
real(8),dimension(kms:kme),intent(inout)   :: t
!- local variables
real(8),dimension(kms:kme) :: a, b, c, d, ntonc, tt
real(8)                    :: x1, x2, tbottom
integer                    :: k
real(8), dimension(:), allocatable    :: aa, bb, cc, dd, xx
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!- calculate newton cooling term
ntonc(:) = 0.0d0
!ntonc(:) = - t(:)/86400.0d0

!- STEP 1: calculate a, b, c, d coefficients
do k = kms+1, kme-1
   x1 =  0.5d0 * (poro(k-1)+poro(k)) /poro(k)
   x2 =  0.5d0 * (poro(k+1)+poro(k)) /poro(k)

   a(k)   =  ( - km(k-1) / ( dz_t(k-1) * dz(k) ) ) * dt *x1
   c(k)   =  ( - km(k)   / ( dz_t(k)   * dz(k) ) ) * dt *x2
   b(k)   =  1.d0 - a(k)  - c(k)
   d(k)   =  t(k)   + dt * ( ft(k)  + ntonc(k) )
end do
!- boundary condition
!- b.c at bottom
!case of tau from bottom
d(kms+1)   =  d(kms+1) - tau_t * (dt/ dz(kms+1))
b(kms+1)   =  1.d0 - c(kms+1)
a(kms+1)   =  0.d0
!- no flux from kme
d(kme-1)     =  d(kme-1) - c(kme-1) *  t(kme)
c(kme-1)     =  0.d0

!- STEP 2: solve simultaneous equations by Thormat method
allocate( aa(kms+1:kme-1))
allocate( bb(kms+1:kme-1))
allocate( cc(kms+1:kme-1))
allocate( dd(kms+1:kme-1))
allocate( xx(kms+1:kme-1))

do k=kms+1, kme-1
   aa(k) = a(k) ; bb(k)=b(k) ; cc(k)=c(k); dd(k)=d(k)
end do

call solv_thomas ( kms+1, kme-1, aa, bb, cc, dd, xx)
do k=kms+1, kme-1
   t(k) = xx(k)
end do

deallocate( aa )
deallocate( bb )
deallocate( cc )
deallocate( dd )
deallocate( xx )

return
end subroutine dyn_moist
!========================================================




!===========================================================
subroutine dyn_momentum (                           &
                       dz, dz_t, poro,              & !-in
                       km, fu, fv, tau_u, tau_v,    & !-in
                       u, v                         & !-inout
                       )
use module_params, only:  kms, kme, dt,sf_surface_physics
implicit none
real(8),dimension(kms:kme),intent(in)     :: dz, dz_t
real(8),dimension(kms:kme),intent(in)     :: poro, km, fu, fv
real(8),dimension(kms:kme),intent(inout)  :: u, v
real(8),intent(in)  :: tau_u, tau_v
!-local variables
real(8),dimension(kms:kme)                :: a, b, c, du, dv
real(8)                                   :: x1, x2
integer                                   :: k
!------------------------------------------------------------
do k = kms+1, kme-1
   x1 =  0.5d0 * (poro(k-1)+poro(k)) /poro(k)
   x2 = 0.5d0 * (poro(k+1)+poro(k)) /poro(k)
   a(k)   =  ( - km(k-1) / ( dz_t(k-1) * dz(k) ) ) * dt * x1
   c(k)   =  ( - km(k)   / ( dz_t(k)   * dz(k) ) ) * dt * x2
   b(k)   =  1.d0 - a(k)  - c(k)
   du(k)  =  u(k) + (fu(k) )*dt
   dv(k)  =  v(k) + (fv(k) )*dt
end do

!- flux from kms
du(kms+1)   =  du(kms+1)  - tau_u * (dt/ dz(kms+1))
dv(kms+1)   =  dv(kms+1)  - tau_v * (dt/ dz(kms+1))
b(kms+1)   =  1.d0 - c(kms+1)
a(kms+1)   =  0.d0

!- no flux from top
du(kme-1)   =  du(kme-1) - c(kme-1) *  u(kme)
dv(kme-1)   =  dv(kme-1) - c(kme-1) *  v(kme)

c(kme-1)   =  0.d0


!- solve tridiag
call solv_tri( a, b, c, du, u)
call solv_tri( a, b, c, dv, v)

return
end subroutine dyn_momentum
!=================================================================








!====================================================
subroutine solv_tri( a, b, c, d, t )
use module_params, only:  kms,kme
implicit none
!   a - sub-diagonal (diagonal below the main diagonal)
!   b - the main diagonal
!   c - sup-diagonal (diagonal above the main diagonal)
!   d - right part
!   t - the answer
!   n = (kme-kms-2) - number of equations
real(8),dimension(kms:kme),intent(in)  :: a,b,c,d
real(8),dimension(kms:kme),intent(out) :: t
real(8),dimension(kms:kme) :: gg,ss
integer :: k


!- initialize g and s
gg(kms+1) = b(kms+1)
ss(kms+1) = d(kms+1)
!- solve for vectors c-prime and d-prime
do k = kms+2,kme-1
   gg(k) = b(k) - a(k)*c(k-1)/gg(k-1)
   ss(k) = d(k) - a(k)*ss(k-1)/gg(k-1)
enddo
!- initialize x
t(kme-1) = ss(kme-1) / gg(kme-1)
!- solve for x from the vectors c-prime and d-prime
do k = kme-2, kms+1, -1
   t(k) = ( ss(k) - c(k) * t(k+1))/gg(k)
end do
return
end subroutine solv_tri
!====================================================




!====================================================
subroutine solv_thomas( ks, ke, a, b, c, d, t )
!====================================================
implicit none
!   a - sub-diagonal (below the main diagonal)
!   b - the main diagonal
!   c - sup-diagonal (above the main diagonal)
!   d - right part
!   t - the answer
!   n = (ke-ks) - number of equations
integer, intent(in)                  :: ks, ke
real(8),dimension(ks:ke),intent(in)  :: a,b,c,d
real(8),dimension(ks:ke),intent(out) :: t
real(8),dimension(ks:ke) :: gg,ss
integer :: k
!----------------------------------------------------

!- initialize g and s
gg(ks) = b(ks)
ss(ks) = d(ks)
!- solve for vectors c-prime and d-prime
do k = ks+1,ke
   gg(k) = b(k) - a(k)*c(k-1)/gg(k-1)
   ss(k) = d(k) - a(k)*ss(k-1)/gg(k-1)
enddo
!- initialize x
t(ke) = ss(ke) / gg(ke)
!- solve for x from the vectors c-prime and d-prime
do k = ke-1, ks, -1
   t(k) = ( ss(k) - c(k) * t(k+1))/gg(k)
end do

return
end subroutine solv_thomas
!=====================================================



end module module_dyn_urban
