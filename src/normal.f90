subroutine normal(pf)
!from physical to normalized ([0-1]) units

use share, only: dp, ndim, n_p, llimits, steps
 
implicit none

!params
real(dp), intent(inout) :: pf(ndim) ! full vector of parameters

!locals
real(dp)                :: ulimit
integer                 :: i

do i=1,ndim
  ulimit=llimits(i)+steps(i)*(n_p(i)-1)
  pf(i)=0.5+(2.*pf(i)-llimits(i)-ulimit)/(2.*(ulimit-llimits(i)))
enddo

end

