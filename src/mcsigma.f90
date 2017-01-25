subroutine	mcsigma(chiscale,w,pf,pf0,obs,lambda_obs,mobs,lsfarr,e,sp,lchi,cov)

!
!	Monte-Carlo error estimates of the error bars.
!
!	We search for solutions after introducing Gaussian errors 
!	in the observations and derive the standard deviation for
!	each parameter.
!

use random
use share, only: dp,ndim,nov,indv,nlambda1,inter,algor,mcruns,mlsf,nlsf
					
implicit none

!in/out
real(dp), intent(in)	   :: chiscale		!chi**2/sum(w_i(f_i-obs_i))^2)
real(dp), intent(in)       :: w(nlambda1)       ! weights
real(dp), intent(in)       :: pf(ndim)          ! vector of fixed parameters
real(dp), intent(in)       :: pf0(ndim)         ! pars read from pfile (in physical units)
real(dp), intent(inout)    :: obs(nlambda1)     ! vector of observations
real(dp), intent(in)       :: lambda_obs(nlambda1)  ! vector of wavelength for observations
real(dp), intent(in)       :: mobs                ! mean or median of obs array
real(dp), intent(in)       :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(in)  	   :: e(nlambda1)		!errors in the observations
real(dp), intent(out)	   :: sp(ndim)		!vector of uncertainties
real(dp), intent(out)  	   :: lchi		!log10(chi**2/(npix-nov+1))
real(dp), intent(out)      :: cov(nov,nov)	!covariance matrix

!locals
real(dp)                   :: pf2(ndim)         !vector of fixed parameters that will get modified 
						!in individual calls to getmin
integer		  	   :: i,j,k		!counters
real(dp)		   :: p(nov)		!vector of var. pars. (1 run)
real(dp)		   :: pm(nov)		!average values
real(dp)		   :: pt(mcruns,nov) !tmp storage for the results
real(dp) 	  	   :: hobs(nlambda1) !holder for the original obs
real(dp)		   :: flux(nlambda1) !interp. flux vector

!initialize
pm(:)=0.0_dp
cov(:,:)=0.0_dp
sp(:)=0.0_dp

!keep the original data; the obs array will be modified here
hobs=obs

!evaluate lchi

call flx(pf,lambda_obs,e,mobs,lsfarr,flux)
lchi=sum(w*(obs(:)-flux(:))**2)
lchi=log10(lchi*chiscale/(nlambda1-nov+1))	

if (lchi >= 4.0d0) then 
!don't even try with lost causes. Very poorly matched spectra 
!are simply assigned sp(*)=-1.00

  	do j=1,nov
		sp(indv(j))=-1.00_dp
  	enddo

else
	
	!initialize the random number generator
	call random_seed()
		
	!repeat mcruns times
	do j=1,mcruns
		
		!introduce Gaussian noise
		do i=1,nlambda1
			obs(i)=random_normal()*e(i)+hobs(i)
		enddo
		
		!initialize
		p(:)=0.5_dp
		!call random_number(p)
		
		!solve
		pf2=pf
		call getmin(algor,w,pf2,pf0,obs,lambda_obs,e,mobs,lsfarr,p)
	
		!store info
		pt(j,:)=p(:)  	!store results to calculate covariance later

		!write(*,*)'mcrun number/sol =',j,p

	enddo
	

	!compute mean
	do j=1,nov
		pm(j)=sum(pt(:,j))/mcruns
	enddo
	!and covariance
	do j=1,nov
	  do k=1,nov
 	   cov(k,j)=sum( (pt(:,k)-pm(k))*(pt(:,j)-pm(j)) )/(mcruns-1._dp)
	  enddo
	enddo

	!use diagonal elements to get std. deviation
  	do j=1,nov
		if (cov(j,j) >=0.0) then
			sp(indv(j))=sqrt(cov(j,j))
		else
			sp(indv(j))=-1._dp
		endif
  	enddo

	
	
endif

!return the original data to obs
obs=hobs

end subroutine	mcsigma
