
subroutine filter2(x,sx,nel,n,y,sy)

! Apply a boxcar filter n+1 pixels wide to a vector x.
! The routine uses linear propagation of errors to derive the errors
! in the smoothed vector from the  errors provided for the original
! data.
!
!
!	IN:	x	Input vector
!		sx	Errors in input vector
!		nel	number of points in the input vectors
!		n	n/2 pixels will be used on each side for even n
!					if not, n/2+1 on the left, 
!					and n/2 on the right
!	
!	OUT:	y	Output vector
!		sy	Resulting vector
!
!

use share, only: dp  

implicit none

!input/output
integer, intent(in)	:: nel,n
real(dp), intent(in)	:: x(nel),sx(nel) !input vector (data and errors)
real(dp), intent(out)	:: y(nel),sy(nel) !output vectors (smooth data and errors)

!locals
integer			:: j,taco,n2,n2t
real(dp)		:: sx2(nel+n)	!=sx**2
real(dp)		:: xx(nel+n),sxx2(nel+n) !tmp copies of x and sx2 with the 
					!extreme values replicated as needed

taco=0
if (n/2. > floor(n/2.)) taco=1
n2=n/2
n2t=n2+taco

do j=1,nel
	sx2(j)=sx(j)*sx(j)
enddo

xx(1:n2t)=x(1)
xx(n2t+1:n2t+nel)=x(1:nel)
xx(n2t+nel+1:n2t+nel+n2)=x(nel)

sxx2(1:n2t)=sx2(1)
sxx2(n2t+1:n2t+nel)=sx2(1:nel)
sxx2(n2t+nel+1:n2t+nel+n2)=sx2(nel)


do j=1,nel
	y(j)=sum(xx(j:j+n2t+n2))/(n+1._dp)
	sy(j)=sqrt(sum(sxx2(j:j+n2t+n2)))/(n+1._dp)
enddo

end subroutine filter2
