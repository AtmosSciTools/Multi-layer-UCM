module module_vars
use module_params, only : kms, kme, nw

implicit none
!- no dimension variables
integer                      :: istep
integer       :: year, month, day, hour, minute, second
!- radiation model variables
real(8)       :: s0                    ! short wave radiation
real(8)       :: elevation             ! solar elevation (radian)
real(8)       :: adi, alpha_s, alpha_c !-solar azimuth: radian, sin, cosin
real(8)       :: sw, sd, ss, lw        !-short (direct, diffuse), long rad
!- surface layer model variables
real(8)       :: tau_u, tau_v, tau_t, tau_qv    ! tau: momentum, heat, moist
real(8)       :: tssl                           ! surface temp. (unused)
!- urban canopy model variables
real(8)       :: hl, wl, rl                     ! urban canopy building index
real(8)       :: svf_sr, svf_rf, svf_wl         !
real(8)       :: svf_r, svf_g
real(8), dimension(:), allocatable        :: svf_w
real(8), dimension(:), allocatable        :: vf_gw
real(8), dimension(:), allocatable        :: vf_wg
real(8), dimension(:,:), allocatable      :: vf_wwp, vf_wwv

real(8)                           :: ts_r, ts_g   ! surfc temp. at roof & road
real(8), dimension(:), allocatable   :: tw_r, tw_g, tgsl  ! inner temp.

real(8), dimension(:,:), allocatable   :: ts   ! surface temp. of walls (2 dim)
real(8), dimension(:,:,:), allocatable :: tw   ! inner temp. of walls ( 4 dim )

integer                                 :: kce
real(8), dimension(:), allocatable      :: z, z_t, dz, dz_t
real(8), dimension(:), allocatable      :: u, v ,t, qv

real(8), dimension(:), allocatable      :: uu, vv ,tt   ! previous step u, v, t
real(8), dimension(:), allocatable      :: p            ! air pressure
real(8), dimension(:), allocatable      :: rho_a, rh    ! air density, rel. hum.
real(8), dimension(:), allocatable      :: km, kh       ! turbulence dif co

real(8), dimension(:), allocatable      :: sw_ra_heating, lw_ra_heating !
real(8), dimension(:), allocatable      :: canop_a, canop_m
real(8), dimension(:), allocatable      :: q, qa, qev   ! building heat

! for input variables
real(8), dimension(:,:), allocatable    :: obs_vars
real(8), dimension(:), allocatable      :: int_vars
character(len=10),allocatable, dimension(:) :: hname
real(8)  :: theta_t, qv_t, u_t, v_t, ws_t
real(8)  :: sw_b, lw_b, sd_b, ss_b

! VARIABLES URBPARM.TBL
real(8)           :: ZR
real(8)           :: ROOF_WIDTH
real(8)           :: ROAD_WITH
real(8)           :: AH
real(8)           :: FRC_URB
real(8)           :: CAPR
real(8)           :: CAPB
real(8)           :: CAPG
real(8)           :: AKSR
real(8)           :: AKSB
real(8)           :: AKSG
real(8)           :: ALBR
real(8)           :: ALBB
real(8)           :: ALBG
real(8)           :: Z0B
real(8)           :: Z0G
real(8)           :: Z0R
real(8)           :: TRLEND
real(8)           :: TBLEND
real(8)           :: TGLEND
real(8)           :: DDZR
real(8)           :: DDZB
real(8)           :: DDZG
real(8), dimension(24)      :: AHDIUPRF
integer           :: THERMAL_INSOL_ROOF
integer           :: THERMAL_INSOL_WALL
integer           :: RLNU
integer           :: BLNU
integer           :: GLNU
integer           :: BOUNDR
integer           :: BOUNDB
integer           :: BOUNDG
integer           :: AHOPTION
integer           :: FRCURBOPT

real(8)    :: COSZ, OMG, DECLIN
real(8)    :: MOLS, MOLR, MOLB, MOLG, MOLC

real(8) :: cpu_sf, cpu_rt, cputime




contains
!===================================================
subroutine var_allocate_memory ( )
allocate( tgsl(1:nw) )
allocate( z(kms:kme+1) )
allocate( dz(kms:kme) )
allocate( z_t(kms:kme) )
allocate( dz_t(kms:kme) )

! physics varibales
allocate( u(kms:kme) )
allocate( v(kms:kme) )
allocate( t(kms:kme) )
allocate( qv(kms:kme) )

allocate( uu(kms:kme) )
allocate( vv(kms:kme) )
allocate( tt(kms:kme) )

allocate( p(kms:kme) )
allocate( rh(kms:kme) )
allocate( rho_a(kms:kme) )

allocate( km(kms:kme) )
allocate( kh(kms:kme) )

allocate( sw_ra_heating(kms:kme) )
allocate( lw_ra_heating(kms:kme) )
allocate( canop_a(kms:kme) )
allocate( canop_m(kms:kme) )
allocate( q(kms:kme) )
allocate( qa(kms:kme) )
allocate( qev(kms:kme) )

return
end subroutine var_allocate_memory









!===================================================
subroutine var_deallocate ( )

deallocate( ts )
deallocate( tw )

deallocate( tw_r )
deallocate( tw_g )
deallocate( tgsl )

deallocate( svf_w )
deallocate( vf_gw )
deallocate( vf_wg )
deallocate( vf_wwp )
deallocate( vf_wwv )

deallocate( z )
deallocate( dz )
deallocate( z_t )
deallocate( dz_t )

deallocate( u )
deallocate( v )
deallocate( t )
deallocate( qv )
deallocate( uu )
deallocate( vv )
deallocate( tt )

deallocate( p )
deallocate( km )
deallocate( kh )
deallocate( rho_a )
deallocate( sw_ra_heating )
deallocate( lw_ra_heating )
deallocate( canop_a )
deallocate( canop_m )
deallocate( q )
deallocate( qa )

return
end subroutine var_deallocate




end module module_vars
