
subroutine getsigma(chiscale,w,pf,obs,lambda_obs,e_obs,mobs,lsfarr,sp,lchi)

!
!
!	Estimate error bars in the fitting parameters
!	The chi**2 projected onto each axis is sampled 
!	with np*2 points between -probe(k) and +probe(k).
!	Returned errors in the parameters correspond to an
!	increase on the chi2 by 1. 
!
!	The optimal parameters enter through the data module share (pf)
!	and are not modified here.
!
	use booklib, only: plot,lsq_fit
	
	use share, only: dp,ndim,nov,indv,npix,nlambda1,probe,badflux,winter,lambda_syn, & 
						mlsf,nlsf
	
	implicit none
	
	!input/output
	real(dp),intent(in)		:: chiscale	!chi**2/sum(w_i(f_i-obs_i))^2)
	real(dp), intent(in)    :: w(nlambda1)           ! weights
	real(dp), intent(in)    :: pf(ndim)              ! vector of fixed parameters
	real(dp), intent(in)    :: obs(nlambda1)         ! vector of observations	
	real(dp), intent(in)    :: lambda_obs(nlambda1)  ! vector of wavelengths for observations
        real(dp), intent(in)    :: e_obs(nlambda1)   ! vector of uncertainties in observations
    real(dp), intent(in)    :: mobs                  ! mean or median of obs array
	real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray
	real(dp), intent(out)	:: sp(ndim) 	!derived uncertainties
								!0=non var. par; -1=unknown
	real(dp), intent(out)	:: lchi		!log10(chi**2) -reduced
	
	!locals
	integer 			:: i,j,k,l
	integer, parameter	:: np=5		!number of points to sample chi2	
	integer, parameter	:: order=2	!order of the polynomial fit
	integer				:: error,nvals
	real(dp)			:: x(np*2),y(np*2),coef(0:order)
	real(dp)			:: chi2,val
	real(dp)			:: flux(nlambda1)
	real(dp)			:: ptmp(ndim)	!temporary storage
	
!	write(*,*)'in getsigma'
    
!	probe contains the size of the (aditive) perturbations we'll 
!	use to observe how chi2 changes as the parameters do


    call flx(pf,lambda_obs,e_obs,mobs,lsfarr,flux)

	lchi=sum(w*(obs(1:nlambda1)-flux(1:nlambda1))**2)
	lchi=log10(lchi*chiscale/(nlambda1-nov+1))	

	sp(:)=0.0_dp		!initialize
	
	if (lchi >= 4.0d0) then 
!	don't even try with lost causes. Very poorly matched spectra 
!	are simply assigned sp(*)=-1.00

	  	do j=1,nov
	  		k=indv(j)
			sp(k)=-1.00_dp
	  	enddo

	else


	  do j=1,nov		!loop over the parameters that were fit

	    	k=indv(j)	!actual dimension in the grid
	    	
	    	y(:)=badflux	!initialize
	    	x(:)=0.0_dp 

				!left side 
	    	do l=1,np	!loop over the grid of points to eval chi**
	    					
	    		
			x(l)=-probe(k)+probe(k)*real(l-1)/real(np)

			!write(*,*)l,x(l)
						
			ptmp=pf
			ptmp(k)=ptmp(k)+x(l)
		
			if (ptmp(k) >= 0.0_dp .and. ptmp(k) <= 1.0_dp) then 
		
				!eval chi**2
				call flx(ptmp,lambda_obs,e_obs,mobs,lsfarr,flux)
				chi2=sum(w*(obs-flux)**2)
				y(l)=chi2*chiscale/(nlambda1-nov+1)
				
			endif
   	   
	    	enddo   	 !end loop over the grid of points to eval chi**2 

				!right side	    	
	    	do l=1,np	!loop over the grid of points to eval chi**
	    					
			x(l+np)=probe(k)*real(l-1)/real(np)
						
			
			ptmp=pf
			ptmp(k)=ptmp(k)+x(l+np)
		
			if (ptmp(k) >= 0.0_dp .and. ptmp(k) <= 1.0_dp) then 
				
				!eval chi**2
				call flx(ptmp,lambda_obs,e_obs,mobs,lsfarr,flux)
				chi2=sum(w*(obs-flux)**2)
				y(l+np)=chi2*chiscale/(nlambda1-nov+1)		
		
			endif
   	   
	    	enddo   	 !end loop over the grid of points to eval chi**2 	    	 
	    	   	   	
	    	!cleanup useless points 
   	   	nvals=count(y > badflux)
   	   	x=pack(x,y > badflux)
   	   	y=pack(y,y > badflux)
   	   	
   	   	if (maxval(y(1:nvals)) == minval(y(1:nvals))) then 
   	   	!if y is flat, then no need to continue, avoiding trouble 
   	   	!with the plotting routine
   	   	
   	   		sp(k)=badflux
   	   		
   	   	else
   	   		
   	   		call plot(y,nvals,6)
			!write(*,*)'plot done'
			!write(*,*)x
			!write(*,*)y
			!write(*,*)nvals,order,error
			!write(*,*)coef
			call lsq_fit(x,y,nvals,order,coef,error)
			!write(*,*)'fit done'
			call plot(coef(0)+coef(1)*x+coef(2)*x**2,nvals,6)
		
			!write(*,*)'coef=',coef
			!write(*,*)'chi**2+1=',(10._dp**lchi+1._dp)
		
			!discriminant

			val=coef(1)**2 - &
			!4._dp*coef(2)*(coef(0)-(10._dp**lchi+change))	
			4._dp*coef(2)*(coef(0)-(10._dp**lchi+1._dp/(nlambda1-nov+1)))

			if (val >= 0) then
				!write(*,*)'(chi**2+1)=',(10._dp**lchi+1._dp)
				!write(*,*)'root1=',(-coef(1)+sqrt(val))/2._dp/coef(2)
				!write(*,*)'root2=',(-coef(1)-sqrt(val))/2._dp/coef(2)	
				sp(k)=	((-coef(1)+sqrt(val))/2._dp/coef(2)-	&
				(-coef(1)-sqrt(val))/2._dp/coef(2))/2._dp
			else
				write(*,*)'negative discriminant - complex roots'
				sp(k)=-1._dp
			endif
		
   	   	   	   	
   	   		!write(*,*)j,pf(k),sp(k)
   	   		
   	   	endif

	   
	  enddo 		!end loop over the parameters	

	endif

	!write(*,*)'exiting getsigma'

	return
	
end subroutine getsigma	
