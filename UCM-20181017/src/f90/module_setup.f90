module module_setup


use module_vars
use module_params
use module_io

contains

subroutine setup ( )

   print*, 'Start module setup', kme
   call read_namelist ( )
   
   call open_files ( )
   call open_nc_file ( )
   call var_allocate_memory ( )
   call print_description( )

   nesting: select case ( ow_nesting )
   case (1)

      write(*,*)  "Nesting with input data"
      write(*,*)  "Input file: ", trim(input_fname)
      call read_obs(input_fname, obs_vars, hname)
      write(*,*) hname, size(hname)

      call check_obs( size(hname), hname,          &
                     sw_opt,lw_opt,sd_opt,ss_opt,  &
                     ws_opt,u_opt,v_opt,           &
                     t_opt,q_opt,ah_opt)
      allocate( int_vars( size(hname) ) )

   case (0) 
      write(*,*)  "No nesting, run ideal mode"     
   case default
      print*, " ERROR :: ow_nesting should be 0 or 1"
      stop
   end select nesting

   print*, '   successed to setup   '
   
return 
end subroutine setup




!=====================================
subroutine read_namelist( )

implicit none
integer                  :: i
integer, parameter       :: im = 6
integer, dimension(im)   :: stat
integer                  :: sumstat 
integer                  :: rstat

print*, 'Reading namelist'
!- open namelist file
open(10, file='./namelist.ucm',action='read',iostat=rstat)
if (rstat /= 0 ) stop 'ERROR: can not open namelist.ucm'

read(unit=10, nml=time_control, iostat=stat(1))
if (stat(1) /= 0 ) print*, 'ERROR: read namelist.ucm time_control'

read(unit=10, nml=domains, iostat=stat(2))
if (stat(2) /= 0 ) print*, 'ERROR: read namelist.ucm domains'

read(unit=10, nml=params, iostat=stat(3)) 
if (stat(3) /= 0 ) print*, 'ERROR: read namelist.ucm params'

read(unit=10, nml=physics, iostat=stat(4)) 
if (stat(4) /= 0 ) print*, 'ERROR: read namelist.ucm physics'

read(unit=10, nml=dynamics, iostat=stat(5))
if (stat(5) /= 0 ) print*, 'ERROR: read namelist.ucm dynamics' 

read (unit=10, NML = nesting, iostat=stat(6) )
if (stat(6) /= 0 ) print*, 'ERROR: read namelist.ucm files.list'


sumstat=0
do i = 1, im 
  sumstat = sumstat + stat(i)
end do
print*, 'Reading namelist'

if ( sumstat /= 0 ) stop  'Error when reading namelist'
close(10)



!==== set run step
nm = int(dfloat((run_day*24+run_hour)*3600)/dt) 
print*, '   run steps: ', nm
history_interval = int(dfloat(output_interval)/dt)
!==================!


return
end subroutine read_namelist
!===============================================================







!================================================================
subroutine print_description( )
implicit none
character(len=100) :: time

write(20,*) '-----------------'
write(20,*) ' physics schemes '
write(20,*) '-----------------'

!write(20,*) " ra_sw_physics :: Kondo94 "
!write(20,*) " ra_lw_physics :: Kondo94 "

pbl_scheme: select case ( bl_pbl_physics )
case (1)
   write(20,*)  " boundary layer physics :: MY_lev2  "
case (0) 
   write(20,*)  " boundary layer physics :: constant  "     
case default
   print*, " ERROR :: Select PBL scheme "
   stop
end select pbl_scheme



suface_scheme: select case ( sf_surface_physics )
case (1)
  write(20,*) " sf_surface_physics :: slab model "  
case (2)
  write(20,*) " sf_surface_physics :: multilayer UCM "

case default
  print*, " ERROR :: Select surface scheme "
  stop
end select suface_scheme

     
return
end subroutine print_description
!=========================================================




!================================================================
subroutine check_obs( nuv, hname,                  &
                     sw_opt,lw_opt,sd_opt,ss_opt,  &
                     ws_opt,u_opt,v_opt,           &
                     t_opt,q_opt,ah_opt)

implicit none
integer, intent(in) :: nuv 
character(len=10), dimension(1:nuv), intent(in) :: hname
integer, intent(inout) :: sw_opt
integer, intent(inout) :: lw_opt
integer, intent(inout) :: sd_opt
integer, intent(inout) :: ss_opt
integer, intent(inout) :: ws_opt
integer, intent(inout) :: u_opt
integer, intent(inout) :: v_opt
integer, intent(inout) :: t_opt
integer, intent(inout) :: q_opt
integer, intent(inout) :: ah_opt
integer :: i
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sw_opt = 0
lw_opt = 0
sd_opt = 0
ss_opt = 0
ws_opt = 0
u_opt  = 0
v_opt  = 0
t_opt  = 0
q_opt  = 0
ah_opt = 0
do i=1, nuv 
   print*, hname(i)
   if (trim(hname(i)) == 'swtop') sw_opt = i
   if (trim(hname(i)) == 'lwtop') lw_opt = i
   if (trim(hname(i)) == 'sdtop') sd_opt = i
   if (trim(hname(i)) == 'sstop') ss_opt = i
   if (trim(hname(i)) == 'wstop') ws_opt = i
   if (trim(hname(i)) == 'utop')  u_opt = i
   if (trim(hname(i)) == 'vtop')  v_opt = i
   if (trim(hname(i)) == 'ttop')  t_opt = i
   if (trim(hname(i)) == 'qtop')  q_opt = i
   if (trim(hname(i)) == 'ahtop') ah_opt = i
end do

!print*, sw_opt,lw_opt,sd_opt,ss_opt,  &
!        ws_opt,u_opt,v_opt,           &
!        t_opt,q_opt,ah_opt

end subroutine check_obs
!================================================================




!================================================================
subroutine read_obs(fname, ovars, hname)
implicit none
character(len=200), intent(in)    :: fname
real(8), dimension(:,:),allocatable, intent(out)  :: ovars
character(len=10),allocatable, dimension(:), intent(out) :: hname
integer :: nrow, ncol, find
integer               :: i, j, iostatus
character(len=200)    :: header
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

find = 999
! find number of rows
open ( find, file=fname, access='SEQUENTIAL',              &
       status='OLD', action='READ', position='REWIND',     &
       iostat=IOSTATUS )
write(*,*) 'Trying to open ',trim(fname)
if(IOSTATUS>0) stop 'ERROR module setup: open obs file'

nrow = 0
do 
   read(find,'()',end=100)
   nrow = nrow + 1
end do
100 close(find)


! find number of rows
open(find, file=fname)
read(find,"(A200)") header
ncol = 0
do i = 1,len_trim(trim(header))
   if (header(i:i)==",") ncol = ncol + 1
end do
close(find)


! allocate variables
allocate(ovars(ncol+1,nrow-1))
allocate(hname(ncol+1))

open(find, file=fname)

do i = 1, nrow 
   if (i==1) then
      read(find,*) (hname(j),j=1,ncol+1)
   else
      read(find,*) (ovars(j,i-1),j=1,ncol+1)
   endif
end do 

close(find)

write(*,*) 'Number of Obs vars: ', ncol
write(*,*) 'Number of Obs rows: ', nrow - 1

!print*, ovars(:,nrow-1)

end subroutine read_obs


end module module_setup
