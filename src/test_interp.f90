
use booklib, only: interp

implicit none

integer, parameter	:: n=3
real			  	:: x(n),y(n),xx(n)
real 				:: yy(n)
integer				:: i,error
real	            :: result

x(1:3)= (/1.,2.,3./)
y(1:3)= (/1.,1.,1./)

xx(1:3)= (/1.,2.,3./)

do i=1,3 
	call interp(x,y,n,xx(i),result,error)
	write(*,*) xx(i),result,error
enddo

end

