! proogram to interpolate observation data
! code by Doan Quang Van
! date: 2013/12/11


program main

implicit none
  integer                                      :: i, nf, stat
  character(len=*),parameter                   :: dir="./preprocess/kugahara/"
  character(len=*),parameter                   :: prefix='interpo_'
  character(len=*),parameter                   :: path="./kugahara/"
  character(len=50),dimension(100)             :: fname
  character(len=100),allocatable, dimension(:) :: flname
  real(8)                                      :: dt  
  namelist /domains/ dt
!  namelist /flist/ nf, flname
!------------------------------------------


  open(10,file='./namelist.pp')
  read(unit=10, nml=domains, iostat=stat) 
     if (stat /= 0 ) print*, 'ERROR: read namelist.input domains'


fname ="unknow"

fname(1) =   "kugahara_0901_sd.txt"
fname(2) =   "kugahara_0901_lw.txt"
fname(3) =   "kugahara_0901_u.txt"
fname(4) =   "kugahara_0901_temp.txt"
fname(5) =   "kugahara_0901_qv.txt"
fname(6) =   "kugahara_0901_ah.txt"
  


do i=1, 100

  if(trim(fname(i))=="unknow") then
     nf = i -1 
     exit
  end if
  
  call data_pp( path, fname(i), dt)

end do



allocate( flname(1:nf))
do i=1, nf
  flname(i) =   trim(dir)//trim(prefix)//trim(fname(i))
  write(*,*) flname(i)
end do




! write files list 

  open(20, file='files.list')
  write(*,*) 'files number:',nf
  !write(20,nml=flist) 

  write(20,*) '&flist'
  write(20,*) 'NF      =',       nf
  do i=1, nf
    write(20,'(a,i2,a,a,a,a)') 'flname(',i,')  =  ', '"',trim(flname(i)),'"'
  end do              
  write(20,*) '/'










!================================================
   contains

   subroutine data_pp( path, fname, dt )

     implicit none
     character(*),intent(in)    :: fname, path 
     real(8), intent(in)        :: dt

     character(len=*),parameter :: prefix='interpo_'
     character(len=50)          :: ffname
     character(len_trim(path)+len_trim(prefix)+len_trim(fname))      :: nfname
     integer                              :: i, k, nrow, nm, nk
     real(8)                              :: tstep 
     real(8)                              :: dtt
     real(8), allocatable,dimension(:)    :: x, x_interp 
!----------------------------------------------------------------------


   ffname=trim(path)//trim(fname)

  ! count number of rows
   open( 10, file=ffname)
   nrow = 0
   do 
       read(10,'()',end=100)
       nrow = nrow + 1
   end do
   100 close(10)
  ! end count number of rows


   nrow = nrow -1
   nk = int(3600.d0/dt)     ! total
   nm = nk*nrow 
   !print *, 'before:',nrow,'after:' ,nm


 !---allocate array
   allocate(x(nrow+1))
   allocate(x_interp(0:nm))


   ! set new file's name
   nfname = trim(path)//trim(prefix)//trim(fname)
   write(*,*) trim(nfname)


  open( 10, file=ffname)
  open( 11, file=nfname)

  
  do i=1, nrow+1
      read(10,*) tstep, x(i)
    !  write(*,*) tstep, x(i)
  enddo




          x_interp(0) = x(1)
          !write(11,*) 0.d0 , x_interp(0)
          dtt = 1.d0/dfloat(nk)
  do i=1, nrow
      do k=1,  nk
          x_interp(i*k) = x(i)  + (x(i+1)-x(i))  *dtt*dfloat(k)
          write(11,*) dfloat(i)-1.d0 + dfloat(k)*dtt, x_interp(k*i)
      enddo
  enddo



  close(10)
  close(11)


  return
  end subroutine data_pp


end program main

