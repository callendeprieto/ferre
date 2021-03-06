
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nea (p,flux)

!nearest-neighbor interpolation

use share, only: dp, ndim,npix,n_p,ntimes,          &
                    f,                              &
	            f_format,f_access,fmtformat
				    
implicit none

!input/output 
real(dp),intent(in)	:: p(ndim)	!ndim vector of wanted pars [0-1]
real(dp),intent(out)	:: flux(npix)	!n vector of output flux
		
!local variables
real(dp)		:: t(ndim)	!input p's run 0-1, t's run 1-n_p(i)

!write(*,*)'in nea'

	
!	get t's from p's. The p's run between 0-1, while t's run 
!	from the min to the max value of the indices of f (1<ti<n_pi) 
t(1:ndim)=p(1:ndim)*(n_p(1:ndim)-1)+1.0_dp

!load wrk
if (f_access == 0) then
  flux(1:npix)=f(1:npix,dot_product(ntimes,nint(t(1:ndim))-1)+1)
else
  !$omp critical
  if (f_format == 0) then
  	read(10,fmtformat,rec=dot_product(ntimes,nint(t(1:ndim))-1)+1) flux(1:npix)
  else
  	read(10,rec=dot_product(ntimes,nint(t(1:ndim))-1)+1) flux(1:npix)
  endif
  !$omp end critical
endif
			
!write(*,*)'exiting nea'

		
end subroutine nea
