program main
implicit none

character(len=200)    :: fname
character(len=10),allocatable, dimension(:) :: hname
real(8), dimension(:,:), allocatable  :: ovars


ovars = -999.d0
fname =   "data/ucm_data/kugahara_topvars__20050901.csv"
! count number of rows

call read_obs(fname, obs_vars, hname)


contains
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

if( IOSTATUS > 0 ) stop '   ERROR: open observation file'

open( find, file=fname)
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



end program
