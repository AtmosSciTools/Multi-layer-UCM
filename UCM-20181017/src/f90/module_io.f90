module module_io

use module_params, only :         &
    kms, kme, nm, run_day,        &
    start_year, start_month,      &
    start_day, start_hour,        &
    start_minute,                 &
    start_second,                 &
    history_interval,             &
    output_interval,              &
    bl_pbl_physics,               &
    sf_surface_physics,           &
    pb, tb, ganma,                &
    lapse_rate, g, dt,            &
    ncid, z_dimid,                &
    zstg_dimid, rec_dimid,        &
    z_varid,                      &
    ofiname

use module_vars, only: kce, THETA_T, QV_T, U_T, V_T, SW_B, LW_B
use netcdf


implicit none
integer :: time_varid
integer :: start(2), count(2)
TYPE ncvar
   integer,dimension(100) :: vid
   character (len = 50),dimension(100) :: vname
   character (len = 50),dimension(100) :: descpt
   character (len = 50),dimension(100) :: units
   character (len = 50),dimension(100) :: dims
   character (len = 50),dimension(100) :: model
END TYPE ncvar
TYPE(ncvar) :: ovar


!=============
 contains

!==================================================
subroutine open_files ( )
implicit none
character (len = 4)    :: run_year

write(run_year,      '(i4)'  ) start_year
open(20,file='.model_desciption_'//run_year//'.txt')  
write(20,*) ''
write(20,*) '      *** MODEL DESCRIPTION ***'
write(20,*) ''

return
end subroutine open_files
!==================================================



!=========================
subroutine close_files ( )
close(20)
return
end subroutine close_files
!=========================



!===========================================================
subroutine read_registry (nuvar)  ! out
! task: to read number of variables registed in Registry.UCM
! then asign readed information to open_netcdf subroutine
! in order to define var name, dimension ...
implicit none
integer               :: IOSTATUS, allocate_status
integer, intent(out)  :: nuvar
character(len=300)    :: string
integer               :: indx
character(len=50)     :: a1, a2, a3, a4, a5, a6

! Open registry file
open ( 100, file='.Registry.UCM',           &
            access='SEQUENTIAL',           &
            status='OLD',                  &
            action='READ',                 &
            position='REWIND',             &
            iostat=IOSTATUS )

if( IOSTATUS > 0 ) print*, '   ERROR: open Registry.UCM'

! read file Registry
! read each line
! ignore meaningless line
! get necessary information from each line:
!    devide line into string
nuvar = 0
READLOOP: do

   read(100,'(A300)', iostat=iostatus) string
   if (iostatus /= 0) exit READLOOP
   if (string(1:1) == "#") cycle READLOOP
   if (trim(string) == "") cycle READLOOP
   read(string,*,iostat=iostatus) a1, a2, a3, a4, a5, a6
   if (iostatus /= 0) stop  "Error: Check Registry file"

   nuvar = nuvar + 1
   ovar%vname(nuvar)   = a1
   ovar%descpt(nuvar)  = a2
   ovar%units(nuvar)   = a3
   ovar%dims(nuvar)    = a4
   ovar%model(nuvar)   = a5

enddo READLOOP
close(100)
return
end subroutine read_registry
!======================================








!====================================
 subroutine open_nc_file ( )
! task: open netcdf file
! in opened file
! there are
implicit none
integer                            :: ii, dim1d(1), wall_dimid
integer, dimension(:), allocatable :: varid
character (len = 50),dimension(:), allocatable   :: varname
character(len=100)                               :: time_units
integer                                          :: nuvar
character (len=200)  :: outfile
character (len=4)    :: run_year
character (len=2)    :: run_day, run_month

write(run_year, '(i4)') start_year
write(run_month,'(i2.2)') start_month
write(run_day,'(i2.2)') start_day

outfile  = trim(ofiname)//"mucmout_"//run_year//"_"//run_month//"_"//run_day//".nc"
print*, "   OUTPUT FILE: ", TRIM(outfile)
! open netcdf file;
! define dimension, variables
call check( nf90_create(trim(outfile), nf90_clobber, ncid) )
call check( nf90_def_dim(ncid, "time",       NF90_UNLIMITED, rec_dimid) )
call check( nf90_def_dim(ncid, "z_dim",      kme,            z_dimid) )
call check( nf90_def_dim(ncid, "z_dim_stag", kme+1,          zstg_dimid))
call check( nf90_def_var(ncid, "z_dim", NF90_REAL, z_dimid,   z_varid) )
call check( nf90_def_var(ncid, "time",  NF90_REAL, rec_dimid, time_varid))
! Get variables information form registry file
call read_registry(nuvar)   ! < get variables information from registry
! however, because not all of registed
! variales should be output regarding to
! used schemes, i.e, if slab scheme was not
! used, therefore the slab-relate variables
! need to be neglected.
! Therefore, schemes information is necessary.
! For example, surface physics scheme:
! When the variable (for output) is determined.
! define variable in netcdf file
! put attributte: descriptions, units ...
surface_scheme: select case(sf_surface_physics)
case(1)
! If slab model was used, then
do ii = 1, nuvar
   IF(trim(ovar%model(ii)).eq."all" .or. &
      trim(ovar%model(ii)).eq."slab".or. &
      trim(ovar%model(ii)).eq."rad") then

      if(trim(ovar%dims(ii)).eq."kn")  then

          call check( nf90_def_var(ncid,                       &
                      trim(ovar%vname(ii)), NF90_FLOAT,        &
                     (/ z_dimid, rec_dimid /), ovar%vid(ii)) )

      elseif(trim(ovar%dims(ii)).eq."n")  then
          call check( nf90_def_var(ncid, trim(ovar%vname(ii)), &
                      NF90_FLOAT, (/rec_dimid/), ovar%vid(ii)) )

      end if
      call check(nf90_put_att(ncid, ovar%vid(ii) ,   &
                 "description",ovar%descpt(ii) ) )
      call check(nf90_put_att(ncid, ovar%vid(ii) ,   &
                 "units",ovar%units(ii) ) )
      call check(nf90_put_att(ncid, ovar%vid(ii) ,   &
                 "dimension",ovar%dims(ii) ) )

   endif
end do

! in the case of urban canopy model
case(2)
call check( nf90_def_dim(ncid, "surface", 6, wall_dimid) )

DO ii = 1, nuvar
   IF(trim(ovar%model(ii)).eq."all" .or. &
      trim(ovar%model(ii)).eq."ucm" .or. &
      trim(ovar%model(ii)).eq."slab".or. &
      trim(ovar%model(ii)).eq."rad") then


      if(trim(ovar%dims(ii)).eq."kn")  then

         call check( nf90_def_var(ncid, trim(ovar%vname(ii)), NF90_FLOAT, &
                    (/ z_dimid, rec_dimid /), ovar%vid(ii)) )

      elseif(trim(ovar%dims(ii)).eq."n")  then
         call check( nf90_def_var(ncid, trim(ovar%vname(ii)),  &
                           NF90_FLOAT, (/rec_dimid/), ovar%vid(ii)) )


      elseif(trim(ovar%dims(ii)).eq."sn")  then
         call check( nf90_def_var(ncid, trim(ovar%vname(ii)),  &
                           NF90_FLOAT, (/wall_dimid,rec_dimid/), ovar%vid(ii)) )


      end if

      call check(nf90_put_att(ncid, ovar%vid(ii) , "description", &
                              ovar%descpt (ii) ) )
      call check(nf90_put_att(ncid, ovar%vid(ii) , "units",       &
                              ovar%units  (ii) ) )
      call check(nf90_put_att(ncid, ovar%vid(ii) , "dimension",   &
                              ovar%dims   (ii) ) )

   END IF
end do


case default
print*, " ERROR :: Select surface scheme "
stop
end select surface_scheme

write (time_units,'(a,i4.4,a,5(i2.2,a))') 'second since ', start_year, '-',  &
                                           start_month, '-', start_day ,     &
                                      ' ', start_hour , ':',start_minute,    &
                                      ':', start_second

call put_global_att( ncid)
call check(nf90_put_att(ncid, time_varid, "units",           time_units     ) )
call check(nf90_put_att(ncid, time_varid, "calendar",        "gregorian"    ) )
call check(nf90_put_att(ncid, time_varid, "timestep",        dt             ) )
call check(nf90_put_att(ncid, time_varid, "output_interval", output_interval) )
! End define mode.
call check( nf90_enddef(ncid) )
return
end subroutine open_nc_file
!=====================================!








!=====================================
subroutine put_global_att( ncid)

implicit none
integer, intent(in)  :: ncid


call check(nf90_put_att(ncid, NF90_GLOBAL, "MUCM", "Developed by Doan Quang Van") )
call check(nf90_put_att(ncid, NF90_GLOBAL, "ra_sw_physics", "Kondo94") )
call check(nf90_put_att(ncid, NF90_GLOBAL, "ra_lw_physics", "Kondo94") )


pbl_scheme: select case ( bl_pbl_physics )
case (1)
   call check(nf90_put_att(ncid, NF90_GLOBAL, "bl_pbl_physics", "MY_lev2") )
case (0)
   call check(nf90_put_att(ncid, NF90_GLOBAL, "bl_pbl_physics", "constant") )
end select pbl_scheme

suface_scheme: select case ( sf_surface_physics )
case (1)
   call check(nf90_put_att(ncid, NF90_GLOBAL, "sf_surface_physics", "slab model") )
case (2)
   call check(nf90_put_att(ncid, NF90_GLOBAL, "sf_surface_physics", "multilayer UCM") )
end select suface_scheme

return
end subroutine put_global_att
!=========================================================










!===================================
subroutine close_nc_file ( )
call check( nf90_close(ncid) )
return
end subroutine close_nc_file
!===================================










!====================================================================
subroutine put_var_nc ( istep, vna, val)

implicit none
integer, intent(in)          :: istep  ! time step
character(len=*),intent(in)  :: vna    ! var name (defined)
real(8),intent(in)           :: val    ! varlues
integer                      :: rec, ii
! NOTE: use global var: history_interval
!                       ovar%name
if( mod(istep,history_interval )== 0) then
    rec = int(istep/history_interval)
    ! Loop for all possible variable then get var ID for favorite var
    do ii = 1, size(ovar%vname)
       if(trim(ovar%vname(ii))==vna) then
          call check( nf90_put_var(ncid, ovar%vid(ii), val, (/rec/)) )
       end if
    end do
end if
return
end subroutine put_var_nc
!=====================================================================




!==============================================================
subroutine put_var2d_nc ( istep, vna, val)
implicit none
integer, intent(in)                :: istep ! time step
character(len=*),intent(in)        :: vna   ! variable name
real(8),dimension(4),intent(in)    :: val   ! variable values
integer  :: rec, ii
! NOTE: use global var: history_interval
!                       ovar%name
if( mod(istep,history_interval )== 0) then
   rec = int(istep/history_interval)
   do ii = 1, size(ovar%vname)
      if(trim(ovar%vname(ii))==vna) then
         call check( nf90_put_var(ncid, ovar%vid(ii), val,  &
                     (/1,rec/),(/6,1/)))
      end if
   end do
end if

return
end subroutine put_var2d_nc
!=============================================================





!======================================
subroutine check(status)
integer, intent ( in) :: status
if(status /= nf90_noerr) then
   print *, trim(nf90_strerror(status))
   stop "Stopped"
end if
end subroutine check
!======================================









!==================================================================
subroutine output(                                    &
                     n,                               &  !-in
                     year, month, day,                &  !-in
                     hour, minute, second,            &  !-in
                     z, z_t, dz, dz_t,                &  !-in
                     canop_a, canop_m,                &  !-in
                     s0, elevation, sw, lw, adi,      &  !-in
                     alpha_s, alpha_c,                &  !-in
                     tau_u, tau_v, tau_t,             &  !-in
                     tau_qv,                          &  !-in
                     u, v, t, qv,                     &  !-in
                     km, kh, q,                       &  !-in
                     uu, vv, tt                       &  !-inout
                  )

implicit none
integer, intent(in)            :: n
integer, intent(in)            :: year, month, day, hour, minute, second
real(8), intent(in)            :: s0,elevation, sw, lw, adi,alpha_s, alpha_c
real(8), intent(in)            :: tau_u, tau_v, tau_t, tau_qv
real(8), dimension(kms:kme+1), intent(in)       :: z
real(8), dimension(kms:kme), intent(in)         :: z_t, dz, dz_t
real(8), dimension(kms:kme), intent(in)         :: canop_a, canop_m
real(8), dimension(kms:kme), intent(in)         :: km, kh, u, v
real(8), dimension(kms:kme), intent(in)         :: t, qv, q
real(8), dimension(kms:kme), intent(inout)      :: uu, vv, tt
real(8)                                         :: tau_tt



!- local variables
integer                                          :: k,i,  dn, ii, rec
real(8)                                          :: hh, dr
real(8), dimension(kms:kme)                      :: theta, spchum, qq
real(8), dimension(kme,10)                      :: ovarval


real(8), dimension(kms:kme)    :: x1, x2, x3, x4, x5
real(8)                        :: sumti,sumto,  sumv, sumu, sumq, sumtau
real(8)                        :: sumqr, sumqw, sumqg
real(8)                        :: xx, xx1, xx2, xx3, xx4, xx5
real(8)                        :: rattg, ratgan , rattr, rattw
real(8)                        :: er, er_rate
integer                        :: ind2m, ind10m
real(8)                        :: THETA2, QV2, U10, V10


!- time by real
hh   =  dfloat(hour)+dfloat(minute)/60.d0 + dfloat(second)/3600.d0
xx1  =  dfloat(n)*dt/86400.d0
dn   =  int(xx1)+1
theta(:)   =  t(:) + tb + ganma*z(:) - 273.15d0  ! or lapse_rate
spchum(:)  =  qv(:) * 1000.d0


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ind2m = 1
! find 2m index and put T2 and Q2
DO K = KMS, KME-1
   IF(Z_T(K) .LE. 2.D0 .AND. Z_T(K+1) .GT. 2.D0) ind2m = K
END DO
!if(N==1) WRITE(*,*) 'INDEX 2M:  ',Z_T(ind2m), Z_T(ind2m+1)
theta2  =  THETA(ind2m)  + ( THETA ( ind2m+1 ) - THETA(ind2m) )  *  &
           (2.D0 - Z_T(ind2m) ) / (Z_T(ind2m+1) - Z_T(ind2m) )

QV2     =  SPCHUM(ind2m) + ( SPCHUM( ind2m+1 ) - SPCHUM( ind2m ) )  &
            * (2.D0 - Z_T( ind2m ) ) / (Z_T(ind2m+1) - Z_T(ind2m) )

! find 10m index and put U10 and V10
ind10m = 1
DO K = KMS, KME-1
   IF(Z_T(K) .LE. 10.D0 .AND. Z_T(K+1) .GT. 10.D0) ind10m = K
enddo
!IF(N==1) WRITE(*,*) 'INDEX 10M:  ',Z_T(ind10m), Z_T(ind10m+1)

U10     =  U(ind10m) + ( U( ind10m+1 ) - U( ind10m ) )  &
           * (2.D0 - Z_T( ind10m ) ) / (Z_T(ind10m+1) - Z_T(ind10m) )

V10     =  V(ind10m) + ( V( ind10m+1 ) - V( ind10m ) )  &
           * (2.D0 - Z_T( ind10m ) ) / (Z_T(ind10m+1) - Z_T(ind10m) )


call put_var_nc ( N, "T2",   THETA2  )
call put_var_nc ( N, "Q2",   QV2     )
call put_var_nc ( N, "U10",  U10     )
call put_var_nc ( N, "V10",  V10     )


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>>>>- write output file every time interval (for 1D variables)
if( n== 1) call check( nf90_put_var(ncid, z_varid, z_t))


count = (/ kme, 1 /)
start = (/ 1, 1 /)
ovarval = 0.d0

if( mod(n,history_interval )== 0) then

   rec = int(n/history_interval)
   start(2) = rec

   call check( nf90_put_var(ncid, time_varid, n*dt, (/rec/)) )

   ovarval(:, 1) = theta(:)
   ovarval(:, 2) = spchum(:)
   ovarval(:, 3) = u(:)
   ovarval(:, 4) = v(:)
   ovarval(:, 5) = km(:)
   ovarval(:, 6) = kh(:)
   ! NOTE that var ID is 1-6 according to Registry.UCM
   do ii = 1, 6
      call check( nf90_put_var(ncid, ovar%vid(ii), &
                  ovarval(:,ii), start = start,    &
                              count = count) )
   end do

end if


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>>>>>>  ---for checking
!- for check heat equation
if(mod(n,history_interval)==0) then
   do k = kms+1 , kme-1
      x3(k)  =  (t(k)-tt(k)) * dz(k) * canop_m(k)     !- [Km]
      x4(k)  =   q(k) * dz(k) * canop_m(k)            !- (Km/s)
      x1(k)  =   ganma * km(k-1) * &
                    (  (canop_m(k+1) + canop_m(k)) - &
                   (canop_m(k) + canop_m(k-1))) * 0.5d0 * dt   !-(Ks/m) !*m2/s
      x2(k)  =   ganma * ( km(k) - km(k-1) )  * dt   !-(Ks/m)    !*m2/s
   end do

   sumto     =  sum(x3(:))                           !- (Km) sum of delta
   sumq      =   sum(x4(:)) * dt                     !- (Km) sum of all q
   sumtau    =  -tau_t * dt  *  canop_m(kms+1)       !- (Km) surface tau
   sumti     =  sumq  +  sumtau  + sum(x2(:))  + sum(x1(:))  !- (Km) sum of input
   er        =  sumto - sumti                                   !- error
   er_rate   =  er / sumto * 100.0d0                            !- error ratio
 end if    !- for check heat equation

! check newton cooling
if(mod(n,8640)==0) then
   do k = kms+1, kme-1
      x1(k) = t(k)*dz(k)*canop_m(k)
   end do
end if

!>>>> - keep record of previous step values of variables
call rec_bfr(                                      &
		n,                                 &  !-in
		u, v, t, km, kh, q, tau_t,         &  !-in
		uu, vv, tt, qq, tau_tt             &  !-out
             )
return
end subroutine output
!==============================================================================







!===============================================================
subroutine rec_bfr(                                 &
		n,                                  &  !-in
		u, v, t, km, kh, q, tau_t,          &  !-in
		uu, vv, tt, qq, tau_tt              &  !-out
                   )

implicit none
integer, intent(in)                          :: n
real(8), intent(in)                          :: tau_t
real(8), dimension(kms:kme), intent(in)      :: km, kh, u, v
real(8), dimension(kms:kme), intent(in)      :: t, q
real(8), dimension(kms:kme), intent(out)     :: uu, vv, tt, qq
real(8)                                      :: tau_tt
integer                                      :: k
!- write files
do k = kms, kme
   uu(k) = u(k)
   vv(k) = v(k)
   tt(k) = t(k)
   qq(k) = q(k)
end do
tau_tt = tau_t
return
end subroutine rec_bfr
!===============================================================

end module module_io
