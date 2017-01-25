subroutine	smooth2(x,sx,nx,n)

! Apply a boxcar filter n+1 pixels wide to the spectra in
! vector x, which may contain several independent wavelength
! ranges + photometry. The error vector is also modified 
! accordingly.
!

use share, only: dp,nsynth,hs
		 		 
implicit none

!input/output
integer, intent(in)	:: n		    !boxcar width is n+1 (see filter2)
integer, intent(in)	:: nx		    !size for x and sx
real(dp),intent(inout) 	:: x(nx),sx(nx) !data and errors

!locals
integer			:: j
real(dp)		:: y(nx),sy(nx) !tmp vars for smoothed data and errs


do j=1,nsynth		
	
	if (abs(hs(j)%res-0._dp) > 1e-6_dp) then	!skip photometry
		
		call filter2(x(hs(j)%pixbegin:hs(j)%pixend),		&
			     sx(hs(j)%pixbegin:hs(j)%pixend),		&
		     		hs(j)%pixend-hs(j)%pixbegin+1,n,	&
		     		y(hs(j)%pixbegin:hs(j)%pixend),	&
		     		sy(hs(j)%pixbegin:hs(j)%pixend))	
		     		
		x(hs(j)%pixbegin:hs(j)%pixend)=y(hs(j)%pixbegin:hs(j)%pixend)
		sx(hs(j)%pixbegin:hs(j)%pixend)=sy(hs(j)%pixbegin:hs(j)%pixend)
		
	endif
enddo

end subroutine smooth2
