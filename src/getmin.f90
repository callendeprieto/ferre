
subroutine getmin(algorithm,irun, opti,fname, chiscale,w,& 
		   pf,pf0,opf,obs,lambda_obs,e_obs,mobs,lsfarr,&
		   incov, spf, lchi, cov)

!
!	Find the minimum the function evaluated in 'objfun'
!	Depending on the value of algorithm, different techniques are used
!	
!  algorithm	method (same codes as shared parameter algor)
!  -1   weighted average over the grid	likely
!   0	pixel with lowest value		minlocus
!   1	Nelder-mead method (Miller's)	minim
!   2   BTR method (Csendes/Miller)	global
!   3   uobyqa method (Powell/Miller)    uobyqa
!   4   truncated Newton  		lmqnbc
!   5   MCMC                             mcmcde



use share, only: dp,ndim,npix,nlambda1,nov,indv,errbar,covprint,     &
			  maxf,iprint,scope,stopcr,nloop,iquad,simp,     &
			  mlsf,nlsf,flen,covprop,minusculo
			  
			  
use nelder_mead,  only	: minim
use btr, only           : global
use uob, only           : uobyqa
use trn, only		: lmqn,lmqnbc
use mcmc, only          : mcmcde
use fun

implicit none

!input/output
integer, intent(in)     :: algorithm    !same codes as 'algor'
integer, intent(in)     :: irun		! index for the run (1...nruns)
integer, intent(in)     :: opti !0=regular run, 1=optimized on
character (len=flen)    :: fname		!spectrum id
real(dp),intent(in)	:: chiscale	!chi**2/sum(w_i(f_i-obs_i))^2)
real(dp), intent(in)    :: w(nlambda1)       ! weights
real(dp), intent(inout) :: pf(ndim)  ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: opf(ndim)     ! pars from pfile (normalized)
real(dp), intent(in)    :: obs(nlambda1)     ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1)         ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1)         ! vector of uncertainties
real(dp), intent(in)    :: mobs              ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf) ! lsfarray
real(dp), intent(in)    :: incov(ndim,ndim)      !input covariance matrix  
						   !used when covprob=1 to propagate covariance in fixed parameters
real(dp), intent(out)	   :: spf(ndim)		!vector of uncertainties
real(dp), intent(out)  	   :: lchi		!log10(chi**2/(nlambda1-nov+1))
real(dp), intent(out)      :: cov(nov,nov)      !covariance matrix


!locals
real(dp) 				:: p(nov)   !vector of variable parameters
real(dp)	   			:: step(nov)	!scope of search for N-M
real(dp)	   			:: var(nov)	!internal errors from N-M
real(dp)				:: flux(nlambda1)
real(dp)	   			:: func		!value of objfun at minimum
integer				:: iter = 0	!actual number of iterations
integer				:: ier = 0  !output error code for minimization routine
integer				:: j,k !dummy indices for loops
real(dp)		  	        :: cof(nov,nov)      !covariance matrix from propagating errors associated to fixed parameters

!specific to BTR
integer, parameter  	:: mglobal = 1		
integer, parameter  	:: n100factor = 100 	!n100=n100factor*nov
integer, parameter  	:: ng00 = 10			!ng0=ng00
integer  				:: nsig					! convergence criterion (number of digits)
integer 				:: n100 !number of sample points to be drawn uniformly in one cycle
integer					:: ng0
integer					:: nc = 0				!number of local minima identified
real(dp)				:: x0(15,20)			!returned solutions
real(dp)				:: f0(20)				!func values at the local minima

!specific to UOBYQA
real(dp)		:: rhobeg      !initial range for parameters (uobyqa)
real(dp)		:: rhoend      !requested final precision (uobyqua)

!specific to TRN
real(dp)		:: g(nov)
integer			:: maxit 
real(dp)		:: eta = 0.25_dp	
real(dp)		:: stepmx = 10._dp
real(dp)		:: faccy, accrcy, xtol
integer			:: ipivot(nov)

!specific to BTR/TRN
real(dp)		:: zeros(nov)  ! min for BTR
real(dp)		:: ones(nov)   ! max for BTR

!specific to MCMC
logical                 :: gr_conv

!set initial values for the search
call pinin(irun,opti,w,pf,opf,obs,lambda_obs,e_obs,mobs,lsfarr,p)


select case (algorithm) 

    case (-1)
		call likely(0,'',w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p)
	        !write(*,*)'back from likely'
        	!write(*,*)'p=',p
     
    case(0) 
        	call minlocus(0,'',w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p)
       	 	!write(*,*)'back from minlocus'
       		!write(*,*)'p=',p
	
    case (1) 
  	    	!write(*,*)'calling minim'
  	    	step(:)=scope
        	call minim(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p, & 
                   step,nov,func,maxf,iprint,stopcr,nloop,iquad,simp,var,ier)
    case (2)
	  	zeros(1:nov)=0.0_dp
	  	ones(1:nov)=1.0_dp
	  	n100=n100factor*nov
	  	ng0=ng00
	  	nc=0
	  	nsig=-log10(stopcr)
	        !$omp critical 
	  	call global(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr, zeros, & 
                   ones,nov,mglobal,n100,ng0,iprint,nsig,x0,nc,f0)
	  	!$omp end critical	  		        
	  	if (nc > 1) then 
	  		write(*,*)' getmin -- multiple minima found, retaining the best only'
	  	endif
	  	p=x0(1:nov,1)
	  	func=f0(1)
    case (3)
	        rhobeg=scope*0.1_dp
	        rhoend=stopcr
                call uobyqa(w,chiscale, pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,nov,p, & 
                   rhobeg,rhoend,iprint,maxf)
                call objfun(w,chiscale, pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,func)
    case (4)
		maxit=nov/2
		eta=0.25_dp
		stepmx=10._dp
		faccy=simp
		accrcy=sqrt(faccy)
		xtol=sqrt(accrcy)
		!call lmqn(w,pf,obs,lambda_obs,e_obs,lsfarr,ier,nov,p,func,g, & 
		!          iprint,maxit,maxf,eta,stepmx,accrcy,xtol)	
		zeros(1:nov)=0.0_dp
	  	ones(1:nov)=1.0_dp
		call lmqnbc(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,ier,nov,p, & 
                   func,g,zeros,ones,ipivot,iprint,maxit,maxf,eta,stepmx,accrcy,xtol)

    case (5)
                call mcmcde(fname,w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr, &
                            p,cov,gr_conv)

		!fill-in error array from cov matrix
		spf(:)=0.0_dp
	  	do j=1,nov
	  		k=indv(j)
			spf(k)=sqrt(cov(j,j))
	  	enddo

		!compute lchi
                call flx(pf,lambda_obs,obs,e_obs,mobs,lsfarr,flux)
		lchi=sum(w*(obs(1:nlambda1)-flux(1:nlambda1))**2)
		lchi=log10(lchi*chiscale/(nlambda1-nov+1))	
	
    case default
        write(*,*) 'ERROR in getmin'
        write(*,*) 'algorithm is not in the allowed range (-1,0,1,2,3,4,5)'
        stop

end select

do j=1,nov
	if (p(j).gt.1.0_dp) then
		pf(indv(j))=1.0_dp
	else if (p(j).le.0.0_dp) then
		pf(indv(j))=1.e-8_dp
	else
		pf(indv(j)) = p(j)
	endif
enddo



!error and cov. matrix determination
if (algorithm /= 5) then 
  lchi=-1.0_dp  
  cov(1:nov,1:nov)=0.0_dp

  select case (errbar)

    case (0)
	call getsigma(chiscale,w,pf,obs,lambda_obs,e_obs,mobs,lsfarr,spf,lchi) 
	if (covprint .eq. 1) call cova(w,lambda_obs,obs,mobs,lsfarr,pf,e_obs,cov)
    case (1)
	call covsigma(chiscale,w,pf,obs,lambda_obs,mobs,lsfarr,e_obs,spf,lchi,cov)	  	
    case (2)
	!we use covsigma to calculate lchi; spf and cov are computed outside getmin
 	call covsigma(chiscale,w,pf,obs,lambda_obs,mobs,lsfarr,e_obs,spf,lchi,cov)
    case (3)
	call pdfsigma(chiscale,w,pf,obs,lambda_obs,mobs,lsfarr,e_obs,spf,lchi,cov)		

    case default
	write(*,*) 'ferre: ERROR'
 	write(*,*) 'errbar=',errbar,' must be 0, 1, 2, or 3'
  	stop
  end select
endif

if (covprop == 1 .and. maxval(incov) > minusculo) then 
        !propagate covariance matrix from fixed parameters
	!inconv=0 indicates this is not needed, and it is used
	!when calling getmin from cofa, to avoid an infinite loop
        call cofa(incov,w,obs,lambda_obs,mobs,lsfarr,pf,e_obs,cof)
        !write(*,*)'cov (original)=',cov
        cov = cov + cof
        !write(*,*)'cof=',cof
        !write(*,*)'cov (after)=',cov

	!use diagonal elements to get std. deviation
  	do j=1,nov
		if (cov(j,j) >=0.0) then
			spf(indv(j))=sqrt(cov(j,j))
		else
			spf(indv(j))=-1._dp
		endif
  	enddo


endif
	
!need to implement a check here on the quality flags from
!the minimiza. routine - if the routine recognizes a failure, this should be
!marked and the solution rejected. 
!some of them don't have such an error flag
       	
      	!write(*,*)'back from minim'
      	!write(*,*)'p=',p
      	!write(*,*)'ier=',ier

end subroutine getmin


