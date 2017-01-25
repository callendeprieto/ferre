subroutine physical(pf)
!from normalized ([0-1]) to physical units

use share, only: dp, ndim, n_p, llimits, steps
 
implicit none

!params
real(dp), intent(inout) :: pf(ndim) ! full vector of parameters

!locals
real(dp)                :: ulimit
integer                 :: i

do i=1,ndim
  ulimit=llimits(i)+steps(i)*(n_p(i)-1)
  pf(i)=pf(i)*(ulimit-llimits(i))+llimits(i)
enddo


end

