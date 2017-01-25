program test_getee

use share, only: inter,ndim,indi,ee
implicit none

integer	:: i

inter=1
ndim=3
indi(1)=3
indi(2)=2
indi(3)=1

call getee

do i=1,(inter+1)**ndim
	write(*,*)i,'->',ee(1,i),ee(2,i),ee(3,i)
enddo

end program test_getee
