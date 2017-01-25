
subroutine filter(x,sx,nel,n,y,sy)

! Apply a boxcar filter n+1 pixels wide to a vector x.
! The routine uses linear propagation of errors to derive the errors
! in the smoothed vector from the  errors provided for the original
! data.
!
!
!	IN:	x	Input vector
!		sx	Errors in input vector
!		nel	number of points in the input vectors
!		n	n/2 pixels will be used on each 
!					if not, n/2 on the left, 
!					and n/2-1 on the right
!	
!	OUT:	y	Output vector
!		sy	Resulting vector
!

use share, only: dp

implicit none


!input/output
integer, intent(in)	:: nel,n
real(dp), intent(in)	:: x(nel),sx(nel) !input vector (data and errors)
real(dp), intent(out)	:: y(nel),sy(nel) !output vectors (smooth data and errors)

!locals
integer			:: j,i1,i2,taco
real(dp)		:: sx2(nel)	  !sx**2

taco=0
if (n/2. > floor(n/2.)) taco=1

do j=1,nel
	sx2(j)=sx(j)*sx(j)
enddo

do j=0,nel-1
	i1=j-n/2
	if (i1 < 0) i1=0
	i2=j+n/2-taco
	if (i2 > nel-1) i2=nel-1
	y(j+1)=sum(x(i1+1:i2+1))/(i2-i1+1._dp)
	sy(j+1)=sqrt(sum(sx2(i1+1:i2+1)))/(i2-i1+1._dp)
enddo

end subroutine filter
