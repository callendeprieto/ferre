MODULE fun

use share, only: dp,ndim,nov,indv,                    &
                 ntie,indtie,typetie,ttie0,ttie,      &
                 nlambda1,mlsf,nlsf,                  &
                 trkout

implicit none

contains

SUBROUTINE objfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, p, func)

 				
implicit none
 
!in/out
real(dp), intent(in)	:: w(nlambda1)		! weights
real(dp), intent(inout)	:: pf(ndim)		! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim)            ! params read from pffile
						! (in physical units)
real(dp), intent(in)	:: obs(nlambda1)	! vector of observations
real(dp), intent(in)	:: lambda_obs(nlambda1)	! vavelengths for observations
real(dp), intent(in)    :: e_obs(nlambda1)      ! uncertainties in observations 
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray
real (dp), intent(in)   :: p(nov)
real (dp), intent(out)  :: func
                                                                                
!locals
real (dp)               :: pp(ndim) 
real (dp)               :: pphys(ndim)          !pp in physical units
integer                 :: i,ii,jj,kk,offset
real(dp)		:: flux(nlambda1)	!interpolated flux vector

pp=pf

!write(*,*)'hola'
do i=1,nov
        pp(indv(i))=p(i)
        !write(*,*)i,p(i)
enddo

if (indtie(1) > 0) then
	pphys(1:ndim)=pp(1:ndim)  
	call physical(pphys)
	if (typetie == 0) then
	  do i=1,ntie
		pphys(indtie(i))=ttie0(i)+sum(ttie(i,1:ndim)*pphys(1:ndim))
	  enddo
	else
	  pphys(1:ndim)=pphys(1:ndim)-pf0(1:ndim)
          do i=1,ntie
                pphys(indtie(i))=ttie0(i)+sum(ttie(i,1:ndim)*pphys(1:ndim))
          enddo
 52       pphys(1:ndim)=pphys(1:ndim)+pf0(1:ndim)
	endif
	call normal(pphys)
	pp(1:ndim)=pphys(1:ndim)
endif

!write(*,*)'pp=',pp
!write(*,*)'obs(1:nlambda1)=',obs(1:nlambda1)
!write(*,*)'flux(1:nlambda1)=',flux(1:nlambda1)
!write(*,*)'w(1:nlambda1)=',w(1:nlambda1)


                                                                                
!evaluate lchi
call flx(pp,lambda_obs,e_obs,mobs,lsfarr,flux)
func=0.0_dp
func=sum(w*(obs(1:nlambda1)-flux(1:nlambda1))**2)
                                                                                
func=func/50._dp

!write(*,*),'func=',func
!write(*,*)'exiting objfun',pp
pf=pp

if (trkout /= 0) then
	if (abs(trkout) == 1) then
		write(11,*) pp
	else 
        	pphys(1:ndim)=pp(1:ndim)
        	call physical(pphys)
		write(11,*) pphys
	endif
	if (trkout < 0) write(12,'(200000(es12.5,1x))') w*(obs-flux)
endif

END SUBROUTINE objfun

!********************************************************************************

function sampler (w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, par_num, zp )


implicit none

!I/O
real(dp), intent(in)	:: w(nlambda1)		! weights
real(dp), intent(inout)	:: pf(ndim)		! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim)            ! params read from pffile
						! (in physical units)
real(dp), intent(in)	:: obs(nlambda1)	! vector of observations
real(dp), intent(in)	:: lambda_obs(nlambda1)	! vavelengths for observations
real(dp), intent(in)    :: e_obs(nlambda1)      ! uncertainties in observations 
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray

integer ( kind = 4 ) par_num
real ( kind = 8 ) sampler
real ( kind = 8 ) zp(par_num)

!locals
integer                 :: i
real (dp)               :: func
real (dp)               :: p(nov)
real (kind = 8)         :: chiscale

p(1:nov)=zp(1:par_num)

call objfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, p, func)

chiscale=sum(w/e_obs**2)/real(nlambda1)

!sampler must be the likelihood = -0.5*chi2
sampler = - 0.5D+00 * func * 50. * chiscale
!sampler = - 0.5D+00 * sum ( ( zp(1:par_num) - 0.5D0 ) ** 2 /0.1d0 **2 )

    
return
end function sampler


END MODULE fun
