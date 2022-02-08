subroutine normal2(pf,ocov)
!from physical to normalized ([0-1]) units

use share, only: dp, ndim, n_p, llimits, steps, ntie, indtie
 
implicit none

!params
real(dp), intent(inout) :: pf(ndim) ! full vector of parameters
real(dp), intent(inout) :: ocov(ndim,ndim)  !full covariance matrix

!locals
real(dp)                :: ulimiti,ulimitj
integer                 :: i,j

do i=1,ndim
  ulimiti=llimits(i)+steps(i)*(n_p(i)-1)
  do j=1,ndim
    ulimitj=llimits(j)+steps(j)*(n_p(j)-1)
    ocov(j,i)=ocov(j,i)/(ulimiti-llimits(i))/(ulimitj-llimits(j))
  enddo
enddo

do i=1,ndim
  ulimiti=llimits(i)+steps(i)*(n_p(i)-1)
  !pf(i)=0.5_dp+(2._dp*pf(i)-llimits(i)-ulimiti)/(2._dp*(ulimit-llimits(i)))
  pf(i)=(pf(i)-llimits(i)) / (ulimiti-llimits(i))
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

end subroutine normal2

