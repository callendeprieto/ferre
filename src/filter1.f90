
subroutine filter1(x,nel,n,y)

! Apply a boxcar filter n+1 pixels wide to a vector x.
!
!
!	IN:	x	Input vector
!		nel	number of points in the input vectors
!		n	n/2 pixels will be used on each side for even n
!					if not, n/2+1 on the left, 
!					and n/2 on the right
!	
!	OUT:	y	Output vector
!

use share, only: dp

implicit none

!input/output
integer, intent(in)	:: nel,n
real(dp), intent(in)	:: x(nel)	  !input vector (data)
real(dp), intent(out)	:: y(nel)	  !output vectors (smoothed data)

!locals
integer			:: j,taco,n2t,n2
real(dp)                :: xx(nel+n) !copy of x with the extremes replicated

taco=0
if (n/2. > floor(n/2.)) taco=1
n2=n/2
n2t=n2+taco

xx(1:n2t)=x(1)
xx(n2t+1:n2t+nel)=x(1:nel)
xx(n2t+nel+1:n2t+nel+n2)=x(nel)
do j=1,nel
	y(j)=sum(xx(j:j+n2t+n2))/(n+1._dp)
enddo

end subroutine filter1
