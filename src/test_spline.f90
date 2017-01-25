
use booklib, only: spline_fit,spline_int

implicit none

integer, parameter	:: n=4,no=11
real			  	:: x(n),y(n),xx(no)
real 				:: yy(n)
integer				:: i,error
real	            :: result
real				:: yp1,ypn,ypp(n)

x(1:4)= (/1.,2.,3.,4./)
y(1:4)= (/1.,2.,3.,4./)

xx(1:no)= (/1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0/)
!xx(1:no)= (/1.0,2.0,3.0/)

call spline_fit(x,y,n,ypp,error)
do i=1,no
	!call interp(x,y,n,xx(i),result,error)
	call spline_int(x,y,n,ypp,xx(i),result,error)
	write(*,*) xx(i),result,error
enddo

end

