program test_getaa

use share, only: aa
implicit none

integer	:: ndim=3,i
integer	:: np(3)

np(1)=2
np(2)=2
np(3)=3

call getaa(ndim,np)

do i=1,2*2*3
	write(*,*)aa(1,i),aa(2,i),aa(3,i)
enddo

end program test_getaa
