module module_dyn_driver

use module_params, only:ganma, tb, ganma, dt,         &
                        kms, kme, sf_surface_physics 

use module_dyn_urban
contains


!==============================================================
subroutine dyn_driver(                                &
                       istep,                         & !-in
                       z, z_t, dz, dz_t,              & !-in
                       canop_a,canop_m, q, qev,       & !-in
                       km, kh,                        & !-in
                       tau_u, tau_v, tau_t, tau_qv,   & !-in
                       u, v, t, qv                    & !-inout
                      )
implicit none
integer,intent(in)                      :: istep
real(8),dimension(kms:kme),intent(in)   :: dz, dz_t, z_t
real(8),dimension(kms:kme+1),intent(in) :: z
real(8),dimension(kms:kme),intent(in)   :: canop_a, canop_m, q, qev
real(8),dimension(kms:kme),intent(in)   :: km, kh
real(8), intent(in)                  :: tau_u, tau_v, tau_t, tau_qv
real(8),dimension(kms:kme),intent(inout):: u, v, t, qv
          
call dyn_urb(                                &
             istep,                          & !-in
             z, z_t, dz, dz_t,               & !-in
             canop_a, canop_m, q, qev,       & !-in
             km, kh,                         & !-in
             tau_u, tau_v, tau_t, tau_qv,    & !-in
             u, v, t, qv                     & !-inout
            )
        
return
end subroutine dyn_driver
!============================================================






!===========================================================
subroutine interpo_data ( istep, vsize, tsize,  & !-i
                          dt_obs, obs_vars, dt, & !-i
                          int_vars )              !-o

implicit none
integer, intent(in)  :: istep, vsize, tsize
real(8), dimension(vsize, tsize), intent(in)  :: obs_vars
real(8), intent(in)  :: dt, dt_obs
real(8), dimension(vsize), intent(out) :: int_vars   
!local variables
integer                      :: n,  i
real(8)                      :: dtt, time
  
!- interpolate observation data 
time =  dfloat(istep)*dt/dt_obs      ! time in real (h)
n    =  int(time) + 1
dtt = time - int(time)

! interpolate for 9 elements
do i=1,  vsize
   int_vars(i) = obs_vars(i,n)   & 
                 + (obs_vars(i,n+1)-obs_vars(i,n))  * dtt        
enddo  
return
end subroutine interpo_data
!=========================================================










!=============================================================
subroutine time_update (                       &
                         istep,                & !- in
                         year, month, day,     & !- inout
                         hour, minute, second  & !- inout
                        )

use module_params, only: dt
implicit none
integer, intent(in)    :: istep
integer, intent(inout) :: year, month, day, hour, minute, second
!-- local variables
integer    :: s1, dm, s2, m1, dh, m2, h1, dd, h2, d1, dmth, d2, &
              mth1, dy, mth2, y1
integer    :: mday1(12), mday2(12)
integer    :: leap

data mday1 / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
data mday2 / 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-leap year distinction
leap = 0
if(mod(year,4).eq.0) leap = 1
if(mod(year,100).eq.0) leap = 0
if(mod(year,400).eq.0) leap = 1
!- new second
s1 = int(dt) + second
dm = int(s1/60)
second = mod(s1,60)
!- new minute
m1 = dm + minute
dh = int(m1/60)
minute = mod(m1,60)
!- new hour
h1 = dh + hour
dd = int(h1/24)
hour = mod(h1,24)

!- new day
d1 = dd + day
dmth = 0
if(leap==0) then       !- case of leap year
   if(d1==mday1(month)+1)  then
      d1   =  1 
      dmth =  1 
   end if
else

   if(d1==mday2(month)+1) then
      d1   =  1 
      dmth =  1 
   end if
end if
!- new month
mth1 = dmth + month
if(mth1==13) then
   mth1= 1
   year  = year + 1
end if


day    = d1
month  = mth1
if(mod(year,4).eq.0) leap = 1
if(mod(year,100).eq.0) leap = 0
if(mod(year,400).eq.0) leap = 1



return
end subroutine time_update
!=========================================================









end module module_dyn_driver
