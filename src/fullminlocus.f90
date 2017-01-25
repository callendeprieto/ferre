subroutine fullminlocus(w,obs,lambda_obs,e_obs,indices)

!
!	When we are searching for all the parameters in a grid
!	simultaneously (ndim=nov), we can take advantage of the intrinsic
!	MINLOC. In this case 'fullminlocus' should be employed instead of
!	minlocus.
!

use share, only:	dp,n_p,ndim,f,ntot,npix,nlambda1,lambda_syn
implicit none

!input/output
real(dp), intent(in)    ::  w(nlambda1)          ! weights
real(dp), intent(in)    ::  obs(nlambda1)        ! vector of observations
real(dp), intent(in)    ::  lambda_obs(nlambda1)        ! vector of wavelengths
real(dp), intent(in)    ::  e_obs(nlambda1)        ! vector of uncertainties
integer, intent(out)	::	indices(ndim)	!minimum location

!locals
real(dp), allocatable   ::      f2(ntot)   		 !temporary storage 
real(dp), allocatable   ::      flux(nlambda1)   !temporary storage 
integer                 ::      j,index(1)


allocate (f2(ntot))	 !allocate f2

do j=1,ntot
	if (winter == 2) then
		call wresample(lambda_syn,f(:,j),npix,lambda_obs,flux,nlambda1)
		f2(j)=sum(w(1:nlambda1)*(flux(1:nlambda1)-obs(1:nlambda1))**2)		
	else
		f2(j)=sum(w(1:nlambda1)*(f(1:nlambda1,j)-obs(1:nlambda1))**2)
	endif
enddo

index=minloc(f2)

end subroutine fullminlocus
