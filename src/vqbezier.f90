
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

	use share, only: dp,npix
	implicit none
	
	integer	::	 i,offset
	real(dp)::	 y1(npix),y2(npix),y3(npix),xp3(3),yp(npix),c0,yprime

	
	if (offset == 0) then

                do concurrent (i=1:npix)
		  yprime=0.5_dp*(y3(i)-y1(i))
		  c0=y2(i)+0.5_dp*yprime
		  yp(i)=y2(i)*xp3(1) + y3(i)*xp3(2) + c0*xp3(3)
                enddo

	else

                do concurrent (i=1:npix)
		  yprime=0.5_dp*(y3(i)-y1(i))
		  c0=y2(i)-0.5_dp*yprime
		  yp(i)=y1(i)*xp3(1) + y2(i)*xp3(2) + c0*xp3(3)
	        enddo	

	endif
		
	end subroutine vqbezier
