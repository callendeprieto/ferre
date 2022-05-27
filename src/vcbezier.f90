
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
        real(dp) :: xp4(4),yp(npix)
        real(dp) :: c0,c1,yprime0,yprime1
        integer  :: i


        do concurrent (i=1:npix)
            yprime0=0.5_dp*(y3(i)-y1(i))
            c0=y2(i)+0.333333333333333333333333_dp*yprime0
            yprime1=0.5_dp*(y4(i)-y2(i))
            c1=y3(i)-0.333333333333333333333333_dp*yprime1
            yp(i)=y2(i)*xp4(1) + y3(i)*xp4(2) + c0*xp4(3) + c1*xp4(4)
        enddo


	end subroutine vcbezier

