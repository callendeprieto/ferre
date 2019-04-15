
subroutine getaa(ndim,np)

! generates a matrix with ndim columns and ntot [=sum_{i=1}^{ndim} np(i)]
! rows which can be used to transform ndim nested loops into a single loop
!
! e.g.
! do i=1,np(1)
!	do j=1,np(2)
!		do k=1,np(3)
!			write(*,*)i,j,k
!		enddo
!	enddo
! enddo
!
! is equivalent to
!
! call getaa(3,np,aa)
! do i=1,np(1)*np(2)*np(3)
!	write(*,*)aa(1,i),aa(2,i),aa(3,i)
! enddo
!
!

use share, only: aa,ntot

implicit none

!input/output
integer, intent(in)			:: ndim
integer, intent(in)			:: np(ndim)


!local
integer					::	i,j
integer					::    	istat !allocate status var
integer					::	v(ndim) !values
integer					::	c(ndim)	!counters
integer					::	nr(ndim)!number of repetitions


if (allocated(aa)) deallocate(aa)
allocate(aa(ndim,ntot),stat=istat)
call checkstat(istat,'aa')

!write(*,*)'ntot=',ntot

!initialize counters and values for the indices
!set number of repetitions 
do i=1,ndim
   v(i)=1
   c(i)=0
   nr(i)=1
   do j=i+1,ndim
   	nr(i)=nr(i)*np(j)
   enddo
   !write(*,*)c(i),v(i),nr(i)
enddo

do i=1,ntot
	do j=1,ndim
		c(j)=c(j)+1
		if (c(j) > nr(j)) then
			v(j)=v(j)+1
			c(j)=1
			if (v(j) > np(j)) then
				v(j)=1
			endif
		endif
		aa(j,i)=v(j)
	enddo
	!write(*,*)aa(:,i)		
enddo

end subroutine getaa
