
	subroutine vcbezier(y1,y2,y3,y4,xp4,yp)

!	Cubic Bezier spline interpolation in 1d
!	y1, y2, y3, and y4 are for x=-1, 0, 1, and 2, respectively
!	the returned value, yp, is the result of the 
!	interpolation for 0<= xp <= 1 
!
! 	The value xp is not actually needed, but 4 derived quantities in xp4
!	xp4=( (1-xp)**3, xp**3, 3*xp*(1-xp)**2, 3*xp**2*(1-xp) )
!
!       This routine is very similar to cbezier, but acts on entire
!       arrays doing the interpolations at various wavelengths concurrently

	use share, only: dp,npix
	implicit none
	
	real(dp) :: y1(npix),y2(npix),y3(npix),y4(npix)
        real(dp) :: xp4(npix,4),yp(npix)
        real(dp) :: c0(npix),c1(npix),yprime0(npix),yprime1(npix)


	yprime0=0.5_dp*(y3-y1)
	c0=y2+0.333333333333333333333333_dp*yprime0
	yprime1=0.5_dp*(y4-y2)
	c1=y3-0.333333333333333333333333_dp*yprime1
		

	!yp=y2*(1.0_dp-xp)**3+y3*xp**3+c0*3.0_dp*xp*(1.0_dp-xp)**2+  &
	!	c1*3.0_dp*xp**2*(1.0_dp-xp)
	
	yp(1:npix)=y2(1:npix)*xp4(1:npix,1) + y3(1:npix)*xp4(1:npix,2) + c0(1:npix)*xp4(1:npix,3) + c1(1:npix)*xp4(1:npix,4)

	end subroutine vcbezier

