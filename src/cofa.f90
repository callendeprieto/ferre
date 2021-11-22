subroutine	cofa(incov,w,obs,lambda_obs,mobs,lsfarr,pf,e,cov)

!
!	Calculation of the covariance matrix (linear term) assuming 
!	Gaussian errors in the observations. See Eq. (14.4.11) of
!	the Numerical Recipes in fortran, and cov=[alpha]**-1. 
!	A solution needs to be determined first and stored in pf, 
!	which enters through the share data module.
!

use booklib, only	: mat_inv
use share, only		: dp,ndim,nov,indv,npix,nlambda1,badflux,winter,lambda_syn, & 
                                    mlsf,nlsf,algor
                                    
implicit none

!in/out
real(dp), intent(in)      :: incov(ndim,ndim) !input covariance matrix
real(dp), intent(in)       :: w(nlambda1)           ! weights
real(dp), intent(in)       :: obs(nlambda1)     ! vector of observations
real(dp), intent(in)       :: lambda_obs(nlambda1)  ! wavelengths for observations
real(dp), intent(in)       :: mobs                 ! mean or median of obs array
real(dp), intent(in)       :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(in)       :: pf(ndim)          ! vector of fixed parameters
real(dp), intent(in)  	    :: e(nlambda1)		!errors in the observations
real(dp), intent(out)      :: cov(nov,nov)	!covariance matrix

!locals
real(dp), parameter	   :: delta=1.1					!factor to compute derivatives
integer		   :: i,j,k						!counters
real(dp)		   :: p(ndim)	        	   		!parameters
real(dp)		   :: flux(nlambda1)				!fluxes
real(dp)		   :: flux1(nlambda1)				!fluxes
real(dp)		   :: flux2(nlambda1)				!fluxes
real(dp)		   :: par1,par2	        			!tmp vars. for partial derivat.
integer		   :: error = 0						!error code from mat_inv
real(dp)		   :: jacob(ndim,ndim)         !jacobian (derivatives of variable wrt fixed quantities)
real(dp) 		   :: pf0(ndim)    	! parameters in physical units
real(dp)		   :: spf(ndim)        ! uncertainties in pf
real(dp)		   :: opf(ndim)        !dummy variable 
real(dp)		   :: lchi             !log of chi**2
real(dp)		   :: fullcov(ndim,ndim) ! full-blown ndim**2 cov. matrix

fullcov(:,:)=0.0

p=pf
opf=pf
pf0=pf
call physical(pf0)

call getmin(algor,1,0,'',1., &
            w,p,pf0,opf,obs,lambda_obs,e,mobs,lsfarr,&
            fullcov,spf,lchi,cov)
	    call flx(p,lambda_obs,e,mobs,lsfarr,flux)

if (flux(1) < 0.0) then 
 cov(:,:)=0.0
else

 do i=1,ndim

   do j=1,ndim


    !only elements in the jacobian that are non-zero are those
    !that correspond to the derivatives of p(i) wrt to p(j) where
    !i corresponds to a variable parameter and j to a fixed one
    if (any(indv == i) .and. .not. any(indv == j)) then 
		
	!two estimates of the partials, increasing
	p(j)=pf(j)*delta
	if (p(j) <= 1.0_dp) then
    		call getmin(algor,1,0,'',1.0, &
               	    w,p,pf0,opf,obs,lambda_obs,e,mobs,lsfarr,&
                	    fullcov,spf,lchi,cov)
	        call flx(p,lambda_obs,e,mobs,lsfarr,flux1)
	endif

	par1=(p(i)-pf(i))/pf(j)/(delta-1._dp)

	!and decreasing the parameters
	p(j)=pf(j)/delta
	if (p(j) >= 0.0_dp) then
    		call getmin(algor,1,0,'',1.0, &
               	    w,p,pf0,opf,obs,lambda_obs,e,mobs,lsfarr,&
                	    fullcov,spf,lchi,cov)
	        call flx(p,lambda_obs,e,mobs,lsfarr,flux2)
	endif
	
	par2=(p(i)-pf(i))/pf(j)/(1._dp/delta-1._dp)

	if (flux1(1) >= 0._dp .and. flux2(1) >=0._dp) then
   		jacob(i,j)=(par1+par2)/2._dp
	else
		jacob(i,j)=0.0_dp
		if (flux1(1) >= 0.0_dp) jacob(i,j)=par1
		if (flux2(1) >= 0.0_dp) jacob(i,j)=par2
	endif


    else
        jacob(i,j)=0.0
    endif

  enddo
	
 enddo

 !write(*,*)'jacob=',jacob
 !write(*,*)'incov=',incov


 fullcov = matmul(jacob,matmul(incov,transpose(jacob)))

 do i=1,nov
   do j=1,nov
     cov(i,j)=fullcov(indv(i),indv(j))
   enddo
 enddo

endif



end subroutine cofa



