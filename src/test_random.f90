
program test_random

use random
implicit none

integer ::	i
real	::	a,b
real	:: c(10)

call random_seed()

do i=1,10
	a=random_normal()
	write(*,*)a
	c(i)=a
enddo

write(*,*)c
write(*,*)c**2

end

