
subroutine fill_ocov (cov,ocov)

!fill-in a ndim x ndim covariance matrix with a nov x nov one
!at the same time normalized parameters are changed to physical

use share, only: dp,ndim,nov,indv,n_p,llimits,steps
			    
implicit none

!input/output 

real(dp),intent(in)	:: cov(nov,nov)
real(dp),intent(out)	:: ocov(ndim,ndim)


!local variables
integer			:: i,l,ii1,ii2		!dummy indices
real(dp)	 	:: ulimit,ulimit2

ocov(:,:)=0.0_dp
do i=1,nov
	ii1=indv(i)
	ulimit=llimits(ii1)+steps(ii1)*(n_p(ii1)-1)
	do l=1,nov
		ii2=indv(l)
		ulimit2=llimits(ii2)+steps(ii2)*(n_p(ii2)-1)
		ocov(ii2,ii1)=cov(l,i)*(ulimit-llimits(ii1))*(ulimit2-llimits(ii2))
	enddo
enddo

		
end subroutine fill_ocov
