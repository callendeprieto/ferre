subroutine nrunsigma (pt,sp,cov)

!this is similar to mcsigma, but here we calculate the cov. matrix from 
!the solutions previously calculated in ferre from abs(nruns)>1
!(no noise has been added to the observations)
!Error bars are calculated from a robust estimate of the 
!standard deviation, and the ratio between these and the plain 
!std. dev. estimates from the sqrt of the diagonal of the 
!covariance matrix is used to scale the whole cov. matrix

use m_mrgrnk
use share, only: dp,ndim,nov,indv,nruns

implicit none

!in/out
real(dp)		   :: pt(nruns,nov) !tmp storage for the results
real(dp), intent(out)	   :: sp(ndim)		!vector of uncertainties
real(dp), intent(out)      :: cov(nov,nov)      !covariance matrix

!locals
integer		  	   :: j,k		!counters
real(dp)		   :: pm(nov)		!average values
integer			   :: index(nruns)	!index for sorting pt(:,j)
real(dp)		   :: scale(nov)	!scaling factors

!initialize
pm(:)=0.0_dp
cov(:,:)=0.0_dp
sp(:)=0.0_dp

!compute mean
do j=1,nov
	pm(j)=sum(pt(:,j))/nruns
enddo

!and covariance
do j=1,nov
 	do k=1,nov
  		cov(k,j)=sum( (pt(:,k)-pm(k))*(pt(:,j)-pm(j)) )/(nruns-1._dp)
	enddo
enddo
!write(*,*)'cov=',cov


!we adopt a robust estimate of the scatter for the error bars
!(which will then become inconsistent with the cov. matrix)
do j=1,nov 
	call mrgrnk(pt(:,j),index)
	sp(indv(j))=(pt(index(int(nruns*0.8415)+1),j)-pt(index(int(nruns*0.1585)+1),j))/2.
enddo
!write(*,*)'robust sp=',sp

!we calculate the scale factors between the original 
!and the robust std. deviation estimates
do j=1,nov
	scale(j)=sp(indv(j))/sqrt(cov(j,j))
enddo

!and use those to scale the whole cov. matrix
do j=1,nov
 	do k=1,nov
  		cov(k,j)=cov(k,j)*scale(k)*scale(j)
	enddo
enddo
!write(*,*)'cov=',cov


!use diagonal elements to get std. deviation
!do j=1,nov
!	if (cov(j,j) >=0.0) then
!		sp(indv(j))=sqrt(cov(j,j))
!	else
!		sp(indv(j))=-1._dp
!	endif
!enddo
!write(*,*)'sp=',sp


end subroutine nrunsigma
