subroutine	smooth1(x,nx,n)

! Apply a boxcar filter n+1 pixels wide to the spectra in
! vector x, which may contain several independent wavelength
! ranges + photometry
!

use share, only: dp,nsynth,hs	
		 
implicit none

!input/output
integer, intent(in)		:: n		  !boxcar width is n+1 (see filter1)
integer, intent(in)		:: nx		  !size of x
real(dp),intent(inout) 	:: x(nx)	  !data



!locals
integer			:: j
real(dp)		:: y(nx)	  !tmp variable for smoothed data



do j=1,nsynth		
	
	if (abs(hs(j)%res-0._dp) > 1e-6_dp) then	!skip photometry
		
		call filter1(x(hs(j)%pixbegin:hs(j)%pixend),		&
		     		hs(j)%pixend-hs(j)%pixbegin+1,n,	&		     		
		     		y(hs(j)%pixbegin:hs(j)%pixend))	

		x(hs(j)%pixbegin:hs(j)%pixend)=y(hs(j)%pixbegin:hs(j)%pixend)		
	endif
enddo


end subroutine smooth1
