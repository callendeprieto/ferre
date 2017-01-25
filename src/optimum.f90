subroutine optimum(obs,lambda_obs,e_obs,mobs,lsfarr,p,w,error)

!
!	Optimizes the weights for the observed spectrum and an
!	an estimate of the model parameters (p)
!
!	On a successful ouput (error=0), the weights have been tuned 
!	and normalized to npix. error=1 signals a failure to optimize
!	and the original weights are not modified.
!
!	Weights are set proportional to the sensitivity of each pixel i,
!	and `sensitivity' is defined as a linear combination of the
!	derivatives of (flux(i)-obs(i))**2 respecto to each of the 
!	parameters p(j). The linear combination may compensate for the
!	overall differences among parameters p(j) in impact on the
!	spectrum. 
!
!	chi(i) is defined as chi(i)=(flux(i)-obs(i))**2
!	X(i,j) is defined as abs(@chi(i)/@p(j))
!	The overall `impact' of a parameter p(j) is defined as
!		Ipct(j)=Sum_i X(i,j) = Sum_i abs(@chi(i)/@p(j))
!	The weights are either set to
!		w(i)=sum(j) X(i,j)/Ipct(j)	if Ipct(j) is not terribly
!	different for all parameters (impact=1), or simply to
!		w(i)=sum(j) X(i,j) (impact=0)
!
!	C. Allende Prieto, July 2004
!

use share, only: dp,ndim,nov,indv,npix,probe,nphotpix,photpixels, & 
			impact,wphot,nlambda1,winter,lambda_syn,              &
			mlsf,nlsf
			
implicit none
	
!input/output
real(dp), intent(in)    ::  obs(nlambda1)    ! vector of observations
real(dp), intent(in)    ::  lambda_obs(nlambda1)  ! vector of wavelengths for observations
real(dp), intent(in)    :: e_obs(nlambda1)        ! vector of uncertainties in observations
real(dp), intent(in)    ::  mobs                 ! mean or median of obs array
real(dp), intent(in)    ::  lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(in)	::	p(ndim)		!estimate of parameters
real(dp), intent(out)   ::  w(nlambda1)      ! weights
integer, intent(out)	::	error		!error flag
	
!locals
integer			::	i,j,k		 			!counters
real(dp)		::	chi,chi2,chi3	 		!chi(i)=(flux(i)-obs(i))**2
real(dp)		:: 	x(nlambda1,nov)	        !X(i,j)=abs(@chi(i)/p(j))
real(dp)		::	Ipct(nov)	 			!vector of `impacts'
real(dp)		::	p2(ndim),p3(ndim)		!aux param vector 
real(dp)		::	flux(nlambda1) 	 !synth flux
real(dp)		::	flux2(nlambda1)	 !aux flux vector
real(dp)		::	flux3(nlambda1)	 !"
real(dp)		::	x2,x3		 !temp storage for elemts. of x
	
error=1

call flx(p,lambda_obs,e_obs,mobs,lsfarr,flux)

!compute x(i,j)
do j=1,nov
	p2=p
	p3=p
	k=indv(j)
	p2(k)=p(k)+probe(k)
	p3(k)=p(k)-probe(k)
	if (p2(k).gt.0.0_dp.and.p2(k).le.1.0_dp.and. &
		p3(k).gt.0.0_dp.and.p3(k).le.1.0_dp) then
		!we have two solutions and take an average
			
			
			call flx(p2,lambda_obs,e_obs,mobs,lsfarr,flux2)
			call flx(p2,lambda_obs,e_obs,mobs,lsfarr,flux2)

			
			do i=1,nlambda1
				chi=(flux(i)-obs(i))**2
				chi2=(flux2(i)-obs(i))**2
				chi3=(flux3(i)-obs(i))**2
				x2=abs(chi2-chi)/probe(k)
				x3=abs(chi3-chi)/probe(k)
				x(i,j)=(x2+x3)/2.0_dp
			enddo
	else
		if (p2(k).gt.0.0_dp.and.p2(k).le.1.0_dp) then
		!only one solution upstream
			
			call flx(p2,lambda_obs,e_obs,mobs,lsfarr,flux2)
			
			do i=1,npix
				chi=(flux(i)-obs(i))**2
				chi2=(flux2(i)-obs(i))**2
				x2=abs(chi2-chi)/probe(k)
				x(i,j)=x2
			enddo		
		
		else
			if (p3(k).gt.0.0_dp.and.p3(k).le.1.0_dp) then
			!only one solution downstream

			call flx(p3,lambda_obs,e_obs,mobs,lsfarr,flux3)

			do i=1,npix
					chi=(flux(i)-obs(i))**2
					chi3=(flux3(i)-obs(i))**2
					x3=abs(chi3-chi)/probe(k)
					x(i,j)=x3
				enddo		
			else
				x(i,j)=0.0_dp
			endif
		endif
	endif
	
enddo	
	

!derive impacts
if (impact.eq.0) then
	Ipct(:)=1.0_dp
else
	do j=1,nov
		Ipct(j)=0.0_dp
		do i=1,npix
			Ipct(j)=Ipct(j)+x(i,j)
		enddo
		write(*,*)j,Ipct(j)
	enddo	
endif


!weights 
do i=1,nlambda1
	x2=0.0_dp
	do j=1,nov
		if (Ipct(j).gt.0.0_dp) x2=x2+x(i,j)/Ipct(j)
	enddo
	if (x2.gt.0.0_dp) then
		w(i)=x2
		error=0	!at least one weight will be updated
	endif
enddo

!reset photometry
do i=1,nphotpix
	w(photpixels(i))=wphot
	write(*,*)'photometry:',photpixels(i),wphot
enddo
	
!and normalize weights 
w=w/sum(w)*real(npix)
	
end subroutine optimum
