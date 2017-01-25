
	subroutine cbezier(y1,y2,y3,y4,xp4,yp)

!	Cubic Bezier spline interpolation in 1d
!	y1, y2, y3, and y4 are for x=-1, 0, 1, and 2, respectively
!	the returned value, yp, is the result of the 
!	interpolation for 0<= xp <= 1 
!
! 	The value xp is not actually needed, but 4 derived quantities in xp4
!	xp4=( (1-xp)**3, xp**3, 3*xp*(1-xp)**2, 3*xp**2*(1-xp) )

	use share, only: dp,mono
	implicit none
	
	real(dp)::	 y1,y2,y3,y4,xp4(4),yp,c0,c1,yprime0,yprime1,mi,ma


	yprime0=0.5_dp*(y3-y1)
	c0=y2+0.333333333333333333333333_dp*yprime0
	yprime1=0.5_dp*(y4-y2)
	c1=y3-0.333333333333333333333333_dp*yprime1
		

        ! force monotonic
	if (mono == 1) then 
	  if (yprime0*yprime1 > 0.0_dp) then
		mi=min(y2,y3)
		ma=max(y2,y3)
		if (c0 < mi) c0=mi
		if (c0 > ma) c0=ma
		if (c1 < mi) c1=mi
		if (c1 > ma) c1=ma		
	  endif
	endif	

	!yp=y2*(1.0_dp-xp)**3+y3*xp**3+c0*3.0_dp*xp*(1.0_dp-xp)**2+  &
	!	c1*3.0_dp*xp**2*(1.0_dp-xp)
	
	yp=y2*xp4(1) + y3*xp4(2) + c0*xp4(3) + c1*xp4(4)

	end subroutine cbezier
