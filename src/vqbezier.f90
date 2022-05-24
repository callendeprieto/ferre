
	subroutine vqbezier(y1,y2,y3,xp3,yp,offset)

!	Quadratic Bezier spline interpolation in 1d
!	y1, y2, and y3 are for x=-1, 0, and 1, respectively
!	the returned value, yp, is the result of the 
!	interpolation for 0<= xp <= 1 if offset=0
!
!	If offset is not 0, it is assumed it is offset=1 and
!	then y1, y2 and y3 are for x=0,1, and 2, respectively
!	and still 0<= xp <= 1
!
!	C. Allende Prieto, Feb 2004
!
! 	The value xp is not actually needed, but 3 derived quantities in xp3
!	xp3=( (1-xp)**2, xp**2, 2*xp*(1-xp) )
!
!       This routine is very similar to qbezier, but acts on entire arrays
! 	doing the interpolations at various wavelengths concurrently

	use share, only: dp, npix
	implicit none
	
	integer	::	offset
	real(dp)	::	y1(npix),y2(npix),y3(npix),xp3(3),yp(npix)
        real(dp)	::      c0(npix),yprime(npix)

	yprime(1:npix)=0.5_dp*(y3(1:npix)-y1(1:npix))
	c0(1:npix)=y2(1:npix)+0.5_dp*yprime(1:npix)
	
	if (offset == 0) then

		!yp=y2*(1.0_dp-xp)**2+y3*xp**2+c0*2.d0*xp*(1.0_dp-xp)
		yp(1:npix)=y2(1:npix)*xp3(1) + y3(1:npix)*xp3(2) + &
                  c0(1:npix)*xp3(3)

	else
		
		!yp=y1*(1.0_dp-xp)**2+y2*xp**2+c0*2.d0*xp*(1.0_dp-xp)
		yp(1:npix)=y1(1:npix)*xp3(1) + y2(1:npix)*xp3(2) + &
                  c0(1:npix)*xp3(3)

	endif
		
	end subroutine vqbezier

