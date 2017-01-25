subroutine lsf_gauss(x,n,fwhm,lsf)

!reconstructs a 1D Gaussian lsf

use share, only: dp, pi

implicit none

!input/output
real(dp), intent(in) :: x(n)    ! pixel location array
integer, intent(in)  :: n       ! number of data points
real(dp), intent(in) :: fwhm    ! fwhm of the lsf
real(dp), intent(out):: lsf(n)  ! lsf profile

!locals
real(dp)     :: sigma_to_fwhm 


sigma_to_fwhm = 2.0_dp*sqrt(-2.0_dp*log(0.5_dp))

lsf(1:n)=1.0_dp/sqrt(2.0_dp*pi)/fwhm*sigma_to_fwhm *  & 
		exp(-(x(1:n)/fwhm*sigma_to_fwhm)**2/2.0_dp)

lsf(1:n)=lsf(1:n)/sum(lsf(1:n))

end subroutine lsf_gauss

