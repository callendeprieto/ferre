subroutine	pdfsigma(chiscale,w,pf,obs,lambda_obs,mobs,lsfarr,e,sp,lchi,cov)

!
!	Calculating the covariance matrix by sampling the likelyhood function
!

use share, only: dp,ndim,nov,indv,npix,nlambda1,lambda_syn,winter, &
                                    mlsf,nlsf


implicit none

!in/out
real(dp), intent(in)	   :: chiscale		!chi**2/sum(w_i(f_i-obs_i))^2)
real(dp), intent(in)       :: w(nlambda1)       ! weights
real(dp), intent(inout)    :: pf(ndim)      ! vector of fixed parameters
real(dp), intent(in)       :: obs(nlambda1)     ! vector of observations
real(dp), intent(in)       :: lambda_obs(nlambda1)  ! wavelengths for observations
real(dp), intent(in)       :: mobs             ! mean or median of obs array
real(dp), intent(in)       :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(in)  	   :: e(nlambda1)		!errors in the observations
real(dp), intent(out)	   :: sp(ndim)		!vector of uncertainties
real(dp), intent(out)  	   :: lchi		!log10(chi**2/(nlambda1-nov+1))
real(dp), intent(out)      :: cov(nov,nov)      !covariance matrix

!locals
integer,parameter  :: nsigma = 4			!how far we go to eval likelyhood
integer, parameter	   :: nsampling	= 1000 			!n. of points for sampling
integer		  	   :: i,j,k		!counters
real(dp)		   :: p(nov)		!vector of var. pars. (1 run)
real(dp)		   :: pm(nov)		!average values
real(dp)		   :: pt(nsampling,nov) !tmp storage for the results
real(dp)		   :: lhood(nsampling)  !likelyhood (exp(-chi2/2))
real(dp)		   :: flux(nlambda1)		!interpolated flux vector
real(dp)		   :: low(nov),high(nov)	!low/hi vals for each norm. parameter
real(dp)		   :: temp						!tmp variable

!evaluate lchi
call flx(pf,lambda_obs,e,mobs,lsfarr,flux)
lchi=0.0_dp
lchi=sum(w*(obs(1:nlambda1)-flux(1:nlambda1))**2)
lchi=log10(lchi*chiscale/(nlambda1-nov+1))	

sp(1:ndim)=0.0_dp		!initialize

!write(*,*)'pdfsigma:  ERROR -- this routine is NOT to be trusted ... yet'
!stop

if (lchi >= 4.0d0) then 
!don't even try with lost causes. Very poorly matched spectra 
!are simply assigned sp(*)=-1.00

  	do j=1,nov
		sp(indv(j))=-1.00_dp
  	enddo

else
	!calculate covariance matrix, just like in covsigma
	call cova(w,lambda_obs,mobs,lsfarr,pf,e,cov)
	
	!use diagonal elements to get std. deviation
  	do j=1,nov
		if (cov(j,j) >=0.0) then
			sp(indv(j))=sqrt(cov(j,j))
		else
			sp(indv(j))=-1._dp
		endif
  	enddo
  	
  	!now we know where to start, let's find the 4-sigma 
  	!intervals in each parameter
  	do j=1,nov
  		low(j)=pf(indv(j))-nsigma*sp(indv(j))
  		high(j)=pf(indv(j))+nsigma*sp(indv(j))
  		if (low(j)   < 0.0) low(j)=0.0_dp
  		if (high(j) >= 1.0) high(j)=0.9999999_dp
  		!write(*,*)'j,low,high=',j,low(j),high(j)
  	enddo

	!generate random sampling data points
	do j=1,nsampling
	   	p(1:ndim)=pf
	   	do i=1,nov
		     	call random_number(temp)
		     	p(indv(i))=temp*(high(i)-low(i))+low(i)
	   	enddo
	   	call flx(p,lambda_obs,e,mobs,lsfarr,flux)
		pt(j,:)=p
	   	temp=0.0_dp
		temp=sum(w*(obs(1:nlambda1)-flux(1:nlambda1))**2)
		temp=temp*chiscale/(nlambda1-nov+1)  !chi**2
		lhood(j)=exp(-temp/2._dp)
	enddo

	!compute weighted mean
	do j=1,nov
		pm(j)=sum(pt(:,j)*lhood)/sum(lhood)
	enddo
	!and weighted covariance
	do j=1,nov
	  do k=1,nov
 	   cov(k,j)=sum( lhood*(pt(:,k)-pm(k))*(pt(:,j)-pm(j)) )/sum(lhood) 
	   ! / 2.3**2
	  enddo
	enddo


!open(30,file='tmp.dat')
!do i=1,nsampling
!	write(30,*) i,chichi(i),lhood(i),pt(i,:)
!	write(*,*),i,nsampling
!enddo
!close(30)

	!update mean values (-- unwanted, the marginalized values are usually biased)
	!pf(1:nov)=pm(1:nov)
	write(*,*)'sum(lhood)=',sum(lhood)
	!write(*,*)'despues  pf=',pf


	!use diagonal elements to get std. deviation
  	do j=1,nov
		if (cov(j,j) >=0.0) then
			sp(indv(j))=sqrt(cov(j,j))
		else
			sp(indv(j))=-1._dp
		endif
  	enddo

endif

end subroutine	pdfsigma
