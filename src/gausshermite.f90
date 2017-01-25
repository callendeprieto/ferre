subroutine gausshermite(x,n,par,npar,ghfun)

!replicates APOGEE's gausshermite.pro

use share, only: dp

!Input/output
real(dp), intent(in)	:: x(n)       ! abscissae
real(dp), intent(in)    :: par(npar)  ! height, center, sigma,H0, H1, H2, H3, H4
integer, intent(in) 	:: n,npar     ! size of x and par
real(dp), intent(out)	:: ghfun(n)   ! Hermite function

!locals
integer, parameter		:: nherm=5
real(dp), parameter		:: pi=3.1415926535897932384626433832795_dp
real(dp), parameter		:: sqr2=sqrt(2._dp),sqr3=sqrt(3._dp)
real(dp), parameter		:: sqr24=sqrt(24._dp),sqrpi=sqrt(pi)
integer					:: i
real(dp)				:: hpar(nherm)
real(dp)				:: hh(nherm),w(n)

hpar(1:nherm)=0.0_dp
if (npar > 3)	hpar(1:npar-3)=par(4:npar)

w(1:n) = (x(1:n) - par(2))/par(3)

hh(1) = hpar(1) - hpar(3)/sqr2 + hpar(5)*3._dp/sqr24
hh(2) = hpar(2) - hpar(4)*3._dp/sqrt(6._dp)
hh(3) = hpar(3)/sqr2 - hpar(5)*(6._dp/sqr24)
hh(4) = hpar(4)/sqrt(6._dp)
hh(5) = hpar(5)/sqr24

ghfun(1:n)=0.0_dp
do i=1,nherm
	ghfun(1:n)= ghfun(1:n) + hh(i)*w(1:n)**(i-1)
enddo
ghfun(1:n)=ghfun(1:n)*exp(-0.5_dp*w(1:n)**2)*par(1)/sqrt(2._dp*pi)

end subroutine gausshermite
