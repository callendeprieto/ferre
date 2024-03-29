subroutine	continuum(x,wx,ox,sx,y,nx,cont,n,rejectcont)

! determine the continuum with one of the following algorithms:
!	cont= 0 -- do nothing
!	cont= 1 -- polynomial fitting
!	cont= 2 -- PEM (split each range in n pieces and divide each subarray by its mean value)
!	cont= 3	-- running mean with a width of n+1 pixels
! 	cont<0  -- do the same as the corresponding positive code but on the ratio x/o 
!

use share, only: dp,nsynth,hs,winter
use booklib
!use lsq
use qr
		 		 
implicit none

!input/output
integer, intent(in)     :: cont		    	!choice of algorithm
integer, intent(in)	:: n		    	!order (cont=1,-1)
					    	!number of pieces per wavelength 
					    	!        segment (cont=2) 
					    	!boxcar width is n+1 (cont=3) 
						!not used for cont=-2 (spline fitting to x/o)
real(dp),intent(in)     :: rejectcont           !rel. error threshold for data rejection when cont=1					    	
integer, intent(in)	:: nx		    	!size for x, wx, ox and sx
real(dp),intent(in) 	:: x(nx),wx(nx),ox(nx), sx(nx)  !data, wavelengths, observations and errors
real(dp),intent(out)    :: y(nx)            	!continuum

!locals
integer			:: p1,p2	    !dummy pixel indices
integer			:: nel              !length of a section
integer			:: nel2             !length of a subsection 
integer			:: j,i		    !loop index
integer			:: error	    !error code for polynomial fit
real(dp)                :: xaxis(nx)        !findgen(nx)
real(dp)		:: w(nx)	    !weights
real(dp),dimension(0:n) :: coef             !coefs. for polynomial fit
real(dp)     		:: xaxis2(nx), x2(nx), ox2(nx), w2(nx)   !temporary arrays for cleaning up 
					    !data  with rel. error > rejectcont
real(dp)                :: ones(nx)         !array of ones to use as weights in polynomial fits
real(dp)		:: ypp(nx)	    !2nd derivatives for cubic splines
real(dp)		:: res		    !temporary variable

ones(:)=1.0_dp

if (cont == 0) then
  y(1:nx)=1._dp
  return
endif


!set abscissae and weights needed for polynomial/spline fitting
if (abs(cont) < 2) then
	w(:)=0._dp
	do i=1,nx
  		xaxis(i)=i*1._dp
		if (sx(i) > 0._dp) w(i)=1._dp/sx(i)**2
	enddo
endif


do j=1,nsynth		
	
	if (abs(hs(j)%res-0._dp) > 1e-6_dp) then	!skip photometry

	  if (winter == 2) then 
		  p1=minval(minloc(wx,wx >= hs(j)%lambdamin))
		  p2=maxval(maxloc(wx,wx <= hs(j)%lambdamax))
	  else
		  p1=hs(j)%pixbegin
		  p2=hs(j)%pixend
	  endif

	  !write(*,*)'lambdamin,lambdamax=',hs(j)%lambdamin,hs(j)%lambdamax
	  !write(*,*)'p1,p2=',p1,p2
	  !write(*,*)'wx(p1),wx(p2)=',wx(p1),wx(p2)

	  nel=p2-p1+1

	  select case (abs(cont))

		  case (1)
		        nel2=0
		        x2(1:nx)=1._dp
		        xaxis2(1:nx)=xaxis(1:nx)
			w2(1:nx)=1._dp
		        ox2(1:nx)=1._dp
		        do i=1,nel
			    if (abs(sx(p1+i-1)/ox(p1+i-1)) < rejectcont) then
			      nel2=nel2+1
			      !xaxis2(nel2)=xaxis(p1+i-1)
			      xaxis2(nel2)=xaxis(i)
			      x2(nel2)=x(p1+i-1)
			      w2(nel2)=w(p1+i-1)
			      ox2(nel2)=ox(p1+i-1)
			    endif
		        enddo
		        !write(*,*) 'rejectcont=',rejectcont, nel, n, nel2
		        if (nel2 <= n) then
		          !give up 
		          write(*,*)'continuum: WARNING'
		          write(*,*) 'Too few points pass the rejectcont=',rejectcont,' filter', nel, n, nel2
		          write(*,*) 'It is being ignored!'
		          nel2=nel
		          xaxis2(1:nel)=xaxis(1:nel)
		          x2(1:nel)=x(p1:p2)
                          w2(1:nel)=w(p1:p2)
		          ox2(1:nel)=ox(p1:p2)
		        endif
		    	if (n == 0) then 
			    !order 0 is just the mean
			    if (cont > 0) then		    
			    	y(p1:p2)=sum(x2(1:nel2))/(nel2*1._dp)
			    else
			    	y(p1:p2)=sum(x2(1:nel2)/ox2(1:nel2))/(nel2*1._dp)
			    endif
			else
			    !fit polynomial
			    coef(:)=0.0_dp
			    if (cont > 0) then 
			    	call polyfit(xaxis2(1:nel2),x2(1:nel2),nel2,n,coef)
			    else
				call polyfit(xaxis2(1:nel2),x2(1:nel2)/ox2(1:nel2),nel2,n,coef)
			    endif

			    !evaluate it
			    y(p1:p2)=coef(0)
			    do i=1,n
				y(p1:p2)=y(p1:p2)+coef(i)*xaxis(1:nel)**i
			    enddo
	  	  	endif

		  case (2)
		    if (cont > 0) then
			call pem(x(p1:p2),nel,n,y(p1:p2))	
		    else
			call pem(x(p1:p2)/ox(p1:p2),nel,n,y(p1:p2))				
		    endif

	  	  case (3)
		    if (cont > 0) then
	    	    	call filter1(x(p1:p2),nel,n,y(p1:p2))
		    else
	    	    	call filter1(x(p1:p2)/ox(p1:p2),nel,n,y(p1:p2))
		    endif

	  	  case default
	    		write(*,*)'ferre: ERROR'
	    		write(*,*)'-- cont must be -3,-2,-1,0,1,2, or 3'
	    		write(*,*)'-- this should have been caught in load_control!'
	    		stop
	  end select
	endif
enddo


end subroutine continuum
