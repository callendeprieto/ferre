subroutine setuplambda(x,nel)

!	Reconstruct a wavelength array for a library or set of
!	libraries from the hs structure.
!

use share, only: dp, long, npix, nsynth, hs

implicit none

!input/output
real(dp), intent(inout) 	:: x(nel)	  !data
integer						:: nel        !number of elements of x (npix or totalnpca)

!locals
integer					:: j
integer(long)           :: i,n

do j=1,nsynth		
	!write(*,*)'setuplambda, j=',j
	!write(*,*)'res=',hs(j)%res
	if (abs(hs(j)%res-0._dp) > 1e-6_dp) then !set wavelengths= 0.0 for photometry
		n=hs(j)%pixend-hs(j)%pixbegin+1
		do i=0,n-1
			x(hs(j)%pixbegin + i) = hs(j)%lambda0 + i*hs(j)%lambda1
		enddo
		if (hs(j)%lws == 1) then 
			x(hs(j)%pixbegin:hs(j)%pixend) =  & 
			10._dp**(x(hs(j)%pixbegin:hs(j)%pixend))
		endif
	else
		x(hs(j)%pixbegin:hs(j)%pixend)= 0._dp
	endif	
enddo

end subroutine setuplambda
