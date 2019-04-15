
subroutine getee

! generates a matrix with ndim columns and base**ndim rows 
! which can be used to transform ndim nested loops (all indices
! running 0 to 1) into a single loop
!
! e.g. for base=2
! do i=0,1
!	do j=0,1
!		do k=0,1
!			write(*,*)i,j,k
!		enddo
!	enddo
! enddo
!
! is equivalent to
! do i=1,2**3
!	write(*,*)ee(1,i),ee(2,i),ee(3,i)
! enddo
!
!
! The matrix eeindi is a copy of ee with its columns
! rearranged according to the interpolation order (indi)
! (they are identical for the default indi values, i.e.
!  for indi = ndim, ndim-1, ndim-2 ...)
!
! The array imap tracks the row in ee which corresponds to
! a given row in eeindi. It allows to sort the wrk array in
! the interpolation routines to match the desired interpolation 
! order.

use share, only: ndim,inter,indi,ee,imap
implicit none

integer			:: i,j,k,nr,c1,c2,base
integer			:: istat	!allocate status var
integer, allocatable 	:: eeindi(:,:)  !copy of ee with columns
					!swapped according to indi
integer 		:: powers(ndim) !(2**0, 2**1, 2**2...)


select case (inter)
	case (1,2,3)
		base=inter+1 !lin,quad., cubic
	case (4)
		base=4       !cubic splines
	case default
		base=2
end select

if (allocated(ee)) deallocate(ee)
allocate (ee(ndim,base**ndim),stat=istat)
call checkstat(istat,'ee')
allocate (eeindi(ndim,base**ndim),stat=istat)
call checkstat(istat,'eeindi')
if (allocated(imap)) deallocate(imap)
allocate (imap(base**ndim),stat=istat)
call checkstat(istat,'imap')

do i=1,ndim
	nr=base**(i-1)
	c1=1
	do j=1,nr
		c2=1
		do k=1,base**ndim/nr
			ee(i,c1)=(c2-1)/base**(ndim-i)
			!write(*,*)i,ee(i,c1)
			c1=c1+1
			c2=c2+1	
		enddo
	enddo
enddo	


!if (indi(1) /= 0) then
	do i=1,ndim
		powers(i)=base**(ndim-i)
	enddo
	
	!copy ee to eeindi and swap columns to match interpolation order
	eeindi=ee
	do i=1,ndim
		eeindi(i,1:base**ndim)=ee(ndim-indi(i)+1,1:base**ndim)
	enddo
	do i=1,base**ndim
		imap(i)=dot_product(eeindi(1:ndim,i),powers)+1
	enddo
!endif


!do i=1,base**ndim
!	write(*,*)'ee=',ee(1:ndim,i)
!enddo
!do i=1,base**ndim
!	write(*,*)'eeindi=',eeindi(1:ndim,i)
!enddo
!
!do i=1,base**ndim
!	write(*,*)'imap=',imap(i)
!enddo

deallocate(eeindi)

end subroutine getee

