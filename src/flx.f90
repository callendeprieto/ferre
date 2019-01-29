
subroutine flx(p,lambda_obs,e_obs,mobs,lsfarr,flux)

!model fluxes derived by 
! a) interpolation 
! b) pca decompression
! c) convolution
! d) resampling in wavelength
! e) square-box filtering
! f) continuum normalization


use share, only: dp,ndim,npix,inter,nfilter,nlambda1,lambda_syn,winter, &
				    npca,totalnpca,pcachi,              &
                                    lsf,mlsf,nlsf,cont,ncont,rejectcont,&
                                    mforce,badflux  
                                    
implicit none

!input/output 
real(dp),intent(in)	:: p(ndim)	!ndim vector of wanted pars [0-1]
real(dp), intent(in)	:: lambda_obs(nlambda1)  ! wavelengths for observations
real(dp), intent(in)	:: e_obs(nlambda1)  ! uncertainties in observations
real(dp), intent(in)    :: mobs ! mean or median of obs array
real(dp), intent(in)	:: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp),intent(out)	:: flux(nlambda1)	!n vector of output flux

!locals
real(dp)		:: sflux(npix)		!synthetic flux vector
real(dp)                :: pflux(totalnpca)     !npca expanded array
real(dp)                :: mflux !mean or median of flux array
real(dp)        	:: cflux(nlambda1)	!continuum


!if outside the 'box', reset fluxes to badflux and quit
if (minval(p).lt.0.0_dp.or.maxval(p).ge.1.0_dp-epsilon(1.0_dp)) then
   flux(1:nlambda1)=badflux
   return
endif	

select case (inter)
case (:0)
	call nea(p,sflux)
case (1)
	call lin(p,sflux)
case (2)
	call qua(p,sflux)
case (3)
	call cub(p,sflux)
case (4)
	call spl(p,sflux)
case default
	 write(*,*) 'flx: ERROR'
	 write(*,*) 'inter=',inter,' must be <= 4'
	 stop
end select
		

if (winter == 2) then 
	if (npca(1) > 0 .and. pcachi == 0) then 
		call decompress(sflux,pflux)
		if (lsf > 0) call convol(pflux,totalnpca,lsfarr,pflux)
		call wresample(lambda_syn,pflux,totalnpca,lambda_obs,flux,nlambda1)
	else
		if (lsf > 0) call convol(sflux,npix,lsfarr,sflux)
		call wresample(lambda_syn,sflux,npix,lambda_obs,flux,nlambda1)
	endif
	if (nfilter >  1) call smooth1(flux,nlambda1,nfilter)
	if (cont > 0) call continuum(flux,lambda_obs,e_obs,cflux,nlambda1, & 
	                             cont,ncont,rejectcont)	
else
	if (npca(1) > 0 .and. pcachi == 0) then !need to eval chi**2 in expanded space
		call decompress(sflux,flux)
		if (lsf > 0) call convol(flux,totalnpca,lsfarr,flux)
		if (nfilter >  1) call smooth1(flux,totalnpca,nfilter)
		if (cont > 0) call continuum(flux,lambda_obs,e_obs,cflux,totalnpca, & 
		                             cont,ncont,rejectcont)
	else
		flux=sflux
		if (lsf > 0) call convol(flux,npix,lsfarr,flux)
		if (nfilter >  1) call smooth1(flux,npix,nfilter)
		if (cont > 0) call continuum(flux,lambda_obs,e_obs,cflux,npix, & 
		                              cont,ncont,rejectcont)
	endif
endif

!write(*,*)'cont=',cont
!write(*,*)'flux=',flux
!write(*,*)'cflux=',cflux

!continuum normalize
if (cont > 0) flux=flux/cflux


!force same mean/median as obs array when mforce>0
if (mforce > 0) then
	if (mforce == 1) then
		mflux=sum(flux(1:nlambda1))/nlambda1
	else 
		call median(flux,nlambda1,mflux)
	endif
	flux=flux-mflux+mobs	
endif

end subroutine flx
