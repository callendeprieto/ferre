subroutine wresample (x,y,n,xx,yy,nn)

use share, only: dp, lambdatol, twinter
use booklib, only: interp, spline_fit, spline_int

implicit none

integer, intent(in)   :: n,nn
real(dp), intent(in)  :: x(n),y(n),xx(nn)
real(dp), intent(out) :: yy(nn)

!locals
integer				  :: i,error
real(dp)              :: result
real(dp)			  :: ypp(n)

select case (twinter)
case (0) !linear wavelength interpolation
  do i=1,nn
	if (xx(i) > lambdatol) then 
		call interp(x,y,n,xx(i),result,error)
		if (error == 0) then
			yy(i)=result
		elseif (error == -1 .and. abs(x(1)-xx(i)) < lambdatol) then
			yy(i)=y(1)
		elseif (error == +1 .and. abs(x(n)-xx(i)) < lambdatol) then
			yy(i)=y(n)
		else
        	write(*,*) 'ERROR in wresample'
        	write(*,*) 'Attempting to interpolate flux to xx=',xx(i)
        	write(*,*) 'but ',x(1),'< x <',x(n)
        	stop
		endif
	else !photometry -- expected at the very beginning of the array
		yy(i)=y(i)
	endif
  enddo
case (1) !cubic splines interpolation
  call spline_fit(x,y,n,ypp,error)
  do i=1,nn
	if (xx(i) > lambdatol) then 
		!write(*,*)'maxval(x)=',maxval(x)
		!write(*,*)'n=',n
		!write(*,*)'xx(i)=',xx(i)
		call spline_int(x,y,n,ypp,xx(i),result,error)
		if (error == 0 .and. xx(i)>x(1)-lambdatol .and. xx(i)<x(n)+lambdatol) then
			yy(i)=result
		else
        	write(*,*) 'ERROR in wresample'
        	write(*,*) 'Attempting to interpolate flux to xx=',xx(i)
        	write(*,*) 'but ',x(1),'< x <',x(n)
        	stop
		endif
	else !photometry -- expected at the very beginning of the array
		yy(i)=y(i)
	endif
  enddo
case default
  write(*,*) 'ERROR in wresample'
  write(*,*) 'twinter must be 0 or 1'
  write(*,*) 'not ',twinter
  stop
end select

end subroutine wresample
