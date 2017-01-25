subroutine likely(wr,fname,w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pos)

!
!	Average estimates
!	pos(i)=sum pos*exp(-chi/min_chi*10) through the whole grid (no interpolation)
!
!	If wr=1 then the chi**2 surface is printed out on fname+'.chi'
!

use share, only: 	dp,flen,liobuffer,n_p,ndim,npix,nlambda1,nov,indv, &
                            llimits,steps,aa,mlsf,nlsf
use fun, only: 	objfun

implicit none

!input/output
integer, intent(in)	::	wr			!write chi**2 file?
character(len=flen),intent(in):: 	fname			   !spectrum id
real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)             ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths for observations
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties for observations
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(out)	:: pos(nov)		   !expected values

!locals
integer			::	i,j,k
integer			::	nt			!number of discrete pts
real(dp)		::	chi
real(dp)		::	min_chi
real(dp)		::	p(nov)			!local vector of pars.
real(dp)		::	ph(nov)			!local vect of phys. par
real(dp)		::	den(nov)		!sum exp(-chi/min_chi*10)
real(dp)		::	ulimit			!temp variable
integer			::	stlen			!tmp counter
character(len=34)	:: 	chifile			!output file 


if (wr > 0) then 
 	stlen=len_trim(fname)
	chifile(1:stlen)=fname
	chifile(stlen+1:stlen+4)='.chi'
	write(*,*)fname
	write(*,*)chifile
	open(5,file=chifile,status='unknown',recl=liobuffer)
endif

!calculate how many discrete points at which to evaluate fitness
nt=n_p(indv(nov))
do j=2,nov
	nt=nt*n_p(indv(nov-j+1))
enddo


!initialize pos and min_chi at center
do j=1,nov
  p(j)=0.5_dp
  pos(j)=p(j)
enddo
call objfun(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,min_chi)

!run over the grid and find the absolute minimum location
do i=1,nt
   do j=1,nov
   	p(j)=(real(aa(j,i))-1._dp)/(real(n_p(indv(j)))-1._dp) !nodes 
   	ulimit=llimits(indv(j))+steps(indv(j))*(n_p(indv(j))-1)
   	ph(j)=p(j)*(ulimit-llimits(indv(j)))+llimits(indv(j)) !physical at nodes
   enddo
   call objfun(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,chi)
   if (wr > 0) write(5,'(1x,21(1x,f8.5))')ph,chi
   if (chi < min_chi) min_chi=chi
enddo

!initialize pos and den 
do j=1,nov
  p(j)=0.0_dp
  pos(j)=0.0_dp
  den(j)=0.0_dp
enddo

!run over the grid and calculate expected values
do i=1,nt
   do j=1,nov
   	p(j)=(real(aa(j,i))-1._dp)/(real(n_p(indv(j)))-1._dp) !nodes 
   enddo
   call objfun(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,chi)
   chi=chi/min_chi
   do j=1,nov
   	pos(j)=pos(j)+p(j)*exp(-chi*10._dp)
   	den(j)=den(j)+exp(-chi*10._dp)
   enddo   
enddo

do j=1,nov
   pos(j)=pos(j)/den(j)
enddo  


if (wr > 0) close(5)

end subroutine likely
