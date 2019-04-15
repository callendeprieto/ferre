program test_lsf_gh

use share, only: dp

implicit none

real(dp)	:: x(100),par(17),lsf(100)
integer		:: i,j

do i=1,100
        x(i)=i-1.0_dp
enddo

open(unit=1,file='kk2')
do i=1,17
read(1,*)par(i)
!write(*,*)par(i)
enddo

call lsf_gh(x,100,50._dp,par,17,lsf)
write(*,*)lsf

end program test_lsf_gh
