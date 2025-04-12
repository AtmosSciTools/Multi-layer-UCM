program main

implicit none

real(8), dimension(1:25)  :: a, b
real(8), dimension(0:96)  :: x
real(8)                   :: h 
integer :: i, k

open (10, file='ah3.2.txt')
open (11, file='kugahara_0901_ah.txt')
      do i=1,25
         read(10,*) a(i), b(i)
         !write(*,'(2f12.6)') a(i), b(i) 
      end do




      do k=0, 3
         h = 1
         do i=1,24
         write(11,'(2f12.6)') float(k)*24. + a(i), b(i) 
       end do
      end do
         write(11,'(2f12.6)') 96., b(25) 



end program main
