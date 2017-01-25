subroutine	cova(w,lambda_obs,mobs,lsfarr,pf,e,cov)

!
!	Calculation of the covariance matrix (linear term) assuming 
!	Gaussian errors in the observations. See Eq. (14.4.11) of
!	the Numerical Recipes in fortran, and cov=[alpha]**-1. 
!	A solution needs to be determined first and stored in pf, 
!	which enters through the share data module.
!

use booklib, only	: mat_inv
use share, only		: dp,ndim,nov,indv,npix,nlambda1,badflux,winter,lambda_syn, & 
                                    mlsf,nlsf
                                    
implicit none

!in/out
real(dp), intent(in)       :: w(nlambda1)           ! weights
real(dp), intent(in)       :: lambda_obs(nlambda1)  ! wavelengths for observations
real(dp), intent(in)       :: mobs                 ! mean or median of obs array
real(dp), intent(in)       :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(in)       :: pf(ndim)          ! vector of fixed parameters
real(dp), intent(in)  	   :: e(nlambda1)		!errors in the observations
real(dp), intent(out)  	   :: cov(nov,nov)	!covariance matrix

!locals
real(dp), parameter	   :: delta=1.1					!factor to compute derivatives
integer		  	   :: i,j,k							!counters
real(dp)		   :: p(ndim)	        	   		!parameters
real(dp)		   :: flux(nlambda1)				!fluxes
real(dp)		   :: flux1(nlambda1)				!fluxes
real(dp)		   :: flux2(nlambda1)				!fluxes
real(dp)		   :: par1,par2	        			!tmp vars. for partial derivat.
real(dp)		   :: partial(nov,nlambda1)				!partial derivatives 
real(dp)		   :: alpha(nov,nov)				!curvature matrix
integer			   :: error = 0						!error code from mat_inv

p=pf
do i=1,nov

	call flx(pf,lambda_obs,e,mobs,lsfarr,flux)
	
	k=indv(i)
	
	!initialize
	flux1(1:nlambda1)=badflux
	flux2(1:nlambda1)=badflux
		
	!two estimates of the partials, increasing
	p(k)=pf(k)*delta
	if (p(k) <= 1.0_dp) then
		call flx(p,lambda_obs,e,mobs,lsfarr,flux1)
	endif

	!and decreasing the parameters
	p(k)=pf(k)/delta
	if (p(k) >= 0.0_dp) then
		call flx(p,lambda_obs,e,mobs,lsfarr,flux2)
	endif
	
		
	do j=1,nlambda1
	
	   par1=(flux1(j)-flux(j))/pf(k)/(delta-1._dp)
	   par2=(flux2(j)-flux(j))/pf(k)/(1._dp/delta-1._dp)

	   if (flux1(j) >= 0._dp .and. flux2(j) >=0._dp) then
	      	partial(i,j)=(par1+par2)/2._dp
	   else
		partial(i,j)=0.0_dp
		if (flux1(j) >= 0.0_dp)	partial(i,j)=par1
		if (flux2(j) >= 0.0_dp) partial(i,j)=par2
		
	   endif
	   
	enddo
enddo

do i=1,nov
	do j=1,nov
		alpha(j,i)=sum(partial(i,1:nlambda1)*partial(j,1:nlambda1)/e(1:nlambda1)**2)
	enddo
enddo

!booklib lib
call mat_inv(alpha,cov,nov,nov,error)
!Pang's lib
!call migs(alpha,nov,cov,flux1)

if (error /= 0) then
	write(*,*) 'cova: ERROR'
	write(*,*) 'mat_inv returned an non-zero error=',error
	write(*,*) 'probably indicating a singular matrix'
	cov(:,:)=-1._dp
endif

end subroutine cova



