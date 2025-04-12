!=================================================================
Program main
!       Fluid Dinamics 1nd dimension simulation
!       Numerics
!          Grid system           : staggard
!          Spatial diffrerence   : 2nd-order centered
!          Time integration      : implicit Euler for diffision
!
!       Physics model
!          PBL                   : Mellor-Yamada level 2
!          Surface layer         : slab model, urban canopy model
!       Code by                  : Doan Quang Van
!                                  Email: doanquangvan@gmail.com
!
!       Date                     : 09 Sep 2012
!       Update 1                 : 15 Feb 2014
!       Update 2                 : 03 Dec 2015
!       Update 3                 : 12 Oct 2018
!==================================================================
use module_vars
use module_setup
use module_initialize
use module_ra_driver
use module_sf_driver
use module_pbl_driver
use module_dyn_driver
use module_urb_ini

implicit none
real(8)                      :: time1, time2, t1_sf, t2_sf

cpu_sf = 0.d0  ! cpu time for surface model
cpu_rt = 0.d0  ! cpu time for ray tracing

call cpu_time(time1)

call setup()

call initialize(                                              &
                 s0, elevation,                               & !- out
                 adi, alpha_s, alpha_c,                       & !- out
                 km, kh,                                      & !- out
                 tau_u, tau_v, tau_t,                         & !- out
                 tau_qv,                                      & !- out
                 z, z_t, dz, dz_t,                            & !- out
                 svf_r, svf_g, svf_w,                         & ! -out
                 vf_gw, vf_wg, vf_wwp, vf_wwv,                & !- out
                 sw, sd, ss, lw,                              & !- out
                 sw_ra_heating, lw_ra_heating,                & !- out
                 kce, canop_a, canop_m, q, qa, qev,           & !- out
                 year, month, day, hour, minute, second,      & !- out
                 tssl, tgsl, ts_g, ts_r, tw_r, tw_g, ts, tw,  & !- out
                 hl, wl, rl,                                  & !- out
                 u, v, p, t,                                  & !- out
                 uu, vv, tt,                                  & !- out
                 qv, rho_a, rh                                & !- out
                )

write(*,*) '   successed to initialize '
write(*,*) '   start time loop '

if(.true.) then

do istep = 1 ,  nm

   ! if run in nesting mode then
   if (ow_nesting==1) then
      call interpo_data ( istep, size(obs_vars,1),  & !-i
                          size(obs_vars,2),         & !-i
                          dt_obs, obs_vars, dt,     & !-i
                          int_vars )

      if (t_opt > 0)  THETA_T =  int_vars(t_opt)
      if (q_opt > 0)  QV_T    =  int_vars(q_opt)
      if (ws_opt > 0) U_T     =  int_vars(ws_opt)
      if (v_opt > 0)  V_T     =  int_vars(v_opt)
      if (u_opt > 0)  V_T     =  int_vars(u_opt)
      if (sw_opt > 0) SW_B    =  int_vars(sw_opt)
      if (lw_opt > 0) LW_B    =  int_vars(lw_opt)
      if (sd_opt > 0) sd_B    =  int_vars(lw_opt)
      if (ss_opt > 0) ss_B    =  int_vars(lw_opt)

      call put_var_nc(ISTEP, "THETA_T", THETA_T )
      call put_var_nc(ISTEP, "QV_T",    QV_T    )
      call put_var_nc(ISTEP, "U_T",     U_T     )
      call put_var_nc(ISTEP, "V_T",     V_T     )
      call put_var_nc(ISTEP, "SW_B",    SW_B    )
      call put_var_nc(ISTEP, "LW_B",    LW_B    )

   endif

   call time_update (                               &
                      istep,                        & !- in
                      year, month, day,             & !- inout
                      hour, minute, second          & !- inout
                     )

   call radiation_driver (                                     &
                       istep, z, dz, dz_t,                     & !- i
                       year, month, day, hour, minute, second, & !- i
                       p, rho_a, t, qv,                        & !- i
                       adi, alpha_s, alpha_c, elevation,       & !- io
                       sw, sd, ss, lw,                         & !- io
                       sw_ra_heating, lw_ra_heating            & !- io
                          )

   if(.true.) then

   call cpu_time(t1_sf)
   call surface_driver (                                      &
                       hour, elevation, alpha_s, alpha_c,     & !- i
                       sw, sd, ss, lw,                        & !- i
                       u, v, p,                               & !- i
                       t, q, qa, qev, tssl, tgsl,             & !- io
                       ts_g, ts_r, tw_r, tw_g, ts, tw,        & !- io
                       qv, rho_a,                             & !- io
                       tau_u, tau_v, tau_t, tau_qv            & !- o
                        )
   call cpu_time(t2_sf)
   cpu_sf = cpu_sf + t2_sf - t1_sf

   call pbl_driver(                       &
	               istep,             &  !-i
                       z, dz, dz_t,       &  !-i
                       u, v, t,           &  !-i
                       km, kh             &  !-o
                    )


   call dyn_driver(                                    &
	               istep,                          & !-i
	               z, z_t, dz, dz_t,               & !-i
	               canop_a, canop_m, q, qev,       & !-i
	               km, kh,                         & !-i
                       tau_u, tau_v, tau_t, tau_qv,    & !-i
                       u, v, t, qv                     & !-io
                   )

   end if
   call output(                                          &
                        istep,                           &  !-in
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

   if(mod(istep,3600*10)==0) write(*,'(a, f10.1)')      &
                       '     time :',dfloat(istep)*dt/3600.


end do ! end time step loop

write(*,'(f10.3)') t(10)

call var_deallocate()
call close_files()
call close_nc_file()
call cpu_time(time2)

write(*,'(a19, f10.5, a5)') "----TOTAL CPU TIME:  ", time2-time1, "SEC"
write(*,'(a19, f10.5, a5)') "----SF MODEL:  ", cpu_sf, "SEC"
write(*,'(a19, f10.5, a5)') "----RT SCHEME:  ", cpu_rt, "SEC"

end if

write(*,*)
write(*,*)'      __________________'
write(*,*)'     |                  | _____'
write(*,*)'     |..... THE END  ...|||"|""|_ '
write(*,*)'     |__________________|||_|___|) '
write(*,*)'     !(@) (@)""""**!(@)(@)****!(@) '
write(*,*)

end program main
