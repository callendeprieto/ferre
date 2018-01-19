
subroutine rmedian(x, y, n, xmed)

! Find the median of a ratio of arrays x/y
! avoiding elements of y with a value of 0.0
!

use share, only: dp

implicit none

INTEGER, INTENT(IN)                    :: n
REAL(dp), INTENT(IN), DIMENSION(n)     :: x,y
REAL(dp), INTENT(OUT)                  :: xmed
REAL(dp), DIMENSION(n)                 :: x2
INTEGER                                :: i,j

i=1
do j=1,n
  if (y(j) /= 0._dp) then
    x2(i)=x(j)/y(j)
    i=i+1
  endif
enddo

call median(x2(1:i-1),i-1,xmed)

end subroutine rmedian
