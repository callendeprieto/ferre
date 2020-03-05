
subroutine chisurf(w,obs,lambda_obs,fname)

!
!	Output the chi**2 surface
!
!	

use share, only: dp,flen,liobuffer,n_p,ntot,ndim,nov,npix,nlambda1,f,winter,lambda_syn
implicit none

!input/output
real(dp), intent(in)        :: w(nlambda1)       ! weights
real(dp), intent(in)        :: obs(nlambda1)     ! vector of observations
real(dp), intent(in)        :: lambda_obs(nlambda1) ! vector of wavelengths for observations
character(len=flen),intent(in):: 	fname		 !spectrum id

!locals
integer			::	i,j,stlen
real(dp)		::	lchi,flux(nlambda1)
character(len=flen)	:: 	chifile		!output file with chi**2 surf

!write(*,*)'in chisurf'


stlen=len_trim(fname)
chifile(1:stlen)=fname(1:stlen)
chifile(stlen+1:stlen+4)='.chi'
write(*,*)fname
write(*,*)chifile
open(5,file=chifile,status='unknown',recl=liobuffer)
write(5,*)fname
write(5,*)obs
do i=1,ntot
	lchi=0.0_dp
	if (winter == 2) then 
		call wresample(lambda_syn,f(j,i),npix,lambda_obs,flux,nlambda1)
		lchi=sum((obs(1:nlambda1)-flux(1:nlambda1))**2)
	else
		lchi=sum((obs(1:nlambda1)-f(1:nlambda1,i))**2)
	endif

	lchi=log10(lchi)
	write(5,*)lchi
enddo
close(5)
!write(*,*)'exiting chisurf'

end subroutine chisurf

