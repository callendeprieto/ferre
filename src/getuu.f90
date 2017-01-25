
subroutine getuu

! generates a matrix (uu, stored in the share module) with nov columns and 
! m=product(indini) rows which can be used to transform nov nested loops into 
! a single loop
!
! e.g. for nov=3 and indini=(2,2,4)
! do i=0,1
!	do j=0,1
!		do k=0,3
!			write(*,*)i,j,k
!		enddo
!	enddo
! enddo
!
! is equivalent to
! do i=1,16
!	write(*,*)uu(1,i),uu(2,i),uu(3,i)
! enddo
!

use share, only: dp,nov,indini,uu

implicit none

integer	    ::	i,j,k,nr,nrd,c1,c2

allocate (uu(nov,product(indini(1:nov))))

do i=1,nov
	c1=1
    nr=product(indini(i:nov))
	!write(*,*)'nr=',nr
	do k=1,product(indini(1:nov))/nr
		c2=1
		nrd=nr/indini(i)
		do j=1,nr
			uu(i,c1)=(c2-1)/nrd
			!write(*,*)'i,c1,uu(i,c1),c2',i,c1,uu(i,c1),c2
			c1=c1+1
			c2=c2+1
		enddo
	enddo
enddo	


end subroutine getuu
