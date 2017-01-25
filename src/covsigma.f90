subroutine	covsigma(chiscale,w,pf,obs,lambda_obs,mobs,lsfarr,e,sp,lchi,cov)

!
!	Calculating the covariance matrix on the assumption of Normal errors
!	and using the diagonal elements to estimate standard errors
!

use share, only: dp,ndim,nov,indv,npix,nlambda1,lambda_syn,winter, &
                                    mlsf,nlsf


implicit none

!in/out
real(dp), intent(in)	   :: chiscale		!chi**2/sum(w_i(f_i-obs_i))^2)
real(dp), intent(in)       :: w(nlambda1)       ! weights
real(dp), intent(in)       :: pf(ndim)      ! vector of fixed parameters
real(dp), intent(in)       :: obs(nlambda1)     ! vector of observations
real(dp), intent(in)       :: lambda_obs(nlambda1)  ! wavelengths for observations
real(dp), intent(in)       :: mobs            ! mean or median of obs array
real(dp), intent(in)       :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(in)  	   :: e(nlambda1)		!errors in the observations
real(dp), intent(out)	   :: sp(ndim)		!vector of uncertainties
real(dp), intent(out)  	   :: lchi		!log10(chi**2/(nlambda1-nov+1))
real(dp), intent(out)      :: cov(nov,nov)      !covariance matrix

!locals
integer		  	   :: i,j						!counters
real(dp)		   :: flux(nlambda1)			!interpolated flux vector

!evaluate lchi
call flx(pf,lambda_obs,e,mobs,lsfarr,flux)
lchi=0.0_dp
lchi=sum(w*(obs(1:nlambda1)-flux(1:nlambda1))**2)
lchi=log10(lchi*chiscale/(nlambda1-nov+1))	

sp(1:ndim)=0.0_dp		!initialize


if (lchi >= 4.0d0) then 
!don't even try with lost causes. Very poorly matched spectra 
!are simply assigned sp(*)=-1.00

  	do j=1,nov
		sp(indv(j))=-1.00_dp
  	enddo

else
	!calculate covariance matrix
	call cova(w,lambda_obs,mobs,lsfarr,pf,e,cov)
	
	!use diagonal elements to get std. deviation
  	do j=1,nov
		if (cov(j,j) >=0.0) then
			sp(indv(j))=sqrt(cov(j,j))
		else
			sp(indv(j))=-1._dp
		endif
  	enddo

	
endif

end subroutine	covsigma
