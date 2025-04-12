module module_params
implicit none

real(8), parameter     ::  g        = 9.8d0 
real(8), parameter     ::  omega    = 7.29d-5        
real(8), parameter     ::  kappa    = 0.4d0
real(8), parameter     ::  sigma    = 5.67d-8
real(8), parameter     ::  cp       = 1004.5d0
real(8), parameter     ::  rd       = 287.0d0   !-
real(8), parameter     ::  lv       = 2.45d6    !-latent heat(J/kg)

!- time_control
integer      :: run_day, run_hour, nm 
integer      :: start_year, start_month, start_day, start_hour,    & 
                start_minute, start_second, history_interval
!- domains
integer            :: kms, kme
real(8)            :: dt
real(8)            :: z_top 
integer            :: grid_type  
!- params
real(8)            :: tb
real(8)            :: ganma 
real(8)            :: lapse_rate        
real(8)            :: lon       
real(8)            :: lat  
real(8)            :: pb       
real(8),parameter  :: humid_urban  = 70.0d0   !-urban humidity(%)
real(8),parameter  :: bdust_urban  = 0.1d0
real(8)            :: Ug
real(8)            :: Vg

!- slab
integer            :: sltype
integer            :: nw
real(8)            :: soil_temp
       
type landuse
   integer          :: lu_id     ! land use index
   real(8)          :: albd      ! aldedo
   real(8)          :: moist     ! moist availability
   real(8)          :: emiss     ! surface emissivity
   real(8)          :: z0        ! surface roughness
   real(8)          :: bo        !
   real(8)          :: lamda
   real(8)          :: cs
   character (len = 100)   :: luname 
end type landuse
type(landuse) ::  lu_slb 
   
 
!- canopy
integer  :: utype 
!- dynamics
integer  :: bl_pbl_physics
integer  :: sf_surface_physics

!- input files
integer :: sw_opt 
integer :: lw_opt 
integer :: sd_opt 
integer :: ss_opt 
integer :: ws_opt 
integer :: u_opt 
integer :: v_opt 
integer :: t_opt 
integer :: q_opt 
integer :: ah_opt 



character(len=100) ::  ofiname
!- Read namelist
integer  :: output_interval

namelist /time_control/  &
     run_day, run_hour, start_year, start_month, start_day,   & 
     start_hour, start_minute, start_second,                  &
     output_interval, ofiname


character(len=200) ::  ics_path
namelist /domains/  &
     dt, kms, kme, z_top, grid_type,ics_path, lon, lat


namelist /params/   & 
     tb, ganma, lapse_rate, pb,  Ug, Vg

namelist /physics/  &
     sf_surface_physics, sltype,  nw , soil_temp, utype

namelist /dynamics/ &
     bl_pbl_physics 

integer   :: ncid
integer   :: z_dimid, z_varid, zstg_dimid, rec_dimid


integer            ::  ow_nesting
real(8)            ::  dt_obs
character(len=200) ::  input_fname
namelist /nesting/  ow_nesting, dt_obs, input_fname


end module module_params
