subroutine normal(pf)
!from physical to normalized ([0-1]) units

use share, only: dp, ndim, n_p, llimits, steps, ntie, indtie
 
implicit none

!params
real(dp), intent(inout) :: pf(ndim) ! full vector of parameters

!locals
real(dp)                :: ulimit
integer                 :: i

do i=1,ndim
  ulimit=llimits(i)+steps(i)*(n_p(i)-1)
  pf(i)=0.5_dp+(2._dp*pf(i)-llimits(i)-ulimit)/(2._dp*(ulimit-llimits(i)))
enddo

if (indtie(1) > 0) then 
  do i=1,ntie
    if (pf(indtie(i)) < 0.0_dp) then
      pf(indtie(i))=0.0_dp
    else
        if (pf(indtie(i)) >= 1.0_dp) pf(indtie(i))=0.999999_dp
    endif
  enddo
endif

end

