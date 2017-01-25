program test_gausshermite

use share, only: dp

implicit none

real(dp) :: x(100),ghfun(100)
real(dp) :: par(8)
integer:: i,j

do i=1,100
	x(i)=real(i)-1.0
enddo

par(1:8)=(/1.0000000  ,     50.000000  ,     1.6529296  ,     1.0000000,   &
-0.019545138 ,    -0.32213661  ,  -0.011828149   ,   0.15485563/)


call gausshermite(x,100,par,8,ghfun)

do i=1,100
write(*,*)ghfun(i)
enddo

end program test_gausshermite
