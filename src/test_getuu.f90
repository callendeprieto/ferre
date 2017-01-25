program test_getuu

use share,only: dp,nov,indini,uu

nov=3

indini(1:nov)=(/2, 1, 8/)
write(*,*)'indini=',indini

call getuu
write(*,*)'back from getuu'

do i=1,product(indini(1:nov))
write(*,*)'uu=',uu(:,i)
enddo

end program test_getuu

