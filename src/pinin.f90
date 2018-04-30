subroutine    pinin(i,opti,w,pf,opf,obs,lambda_obs,e_obs,mobs,lsfarr,p)

!
!  Routine that sets the initial values for
!  the variable parameters in the search. This is done
!  following the directives in the array indini, which
!   can be overriden using the variable init.
!

use share, only: dp,ndim,nov,indv,indini,uu,   &
                 nlambda1,mlsf,nlsf,           &
                 algor,init


implicit none

!in/out
integer, intent(in)   :: i		! index for the run (1...nruns)
integer, intent(in)   :: opti !0=regular run, 1=optimized on
																														!and therefore a previous regular run
																														!has taken place and the results are in p
real(dp), intent(in)  :: w(nlambda1)   ! weights
real(dp), intent(in)  :: pf(ndim)      ! vector of fixed parameters
real(dp), intent(in)  :: opf(ndim)     ! vector of fixed parameters
																																				! before any previous optimization
real(dp), intent(in)  :: obs(nlambda1)	! vector of observations
real(dp), intent(in)  :: lambda_obs(nlambda1)  ! wavelengths for observations
real(dp), intent(in)  :: e_obs(nlambda1)  ! uncertainties for observations
real(dp), intent(in)  :: mobs             ! mean or median of obs array
real(dp), intent(in)  :: lsfarr(mlsf,nlsf)    ! lsfarray
real (dp), intent(inout) :: p(nov)

!locals
integer   :: j,ii,kk   		


if (opti == 0) then !regular run

 !initialization for the amoeba
 select case (init)
	  case (0) !use input values from pfile
		p(1:nov)=pf(indv(1:nov))
	  case (1) !start as coded in indini					
   	 	do j=1,nov
			if (indini(j) == 0) then
				call random_number(p(j))
			else
          		p(j)=uu(j,i)/dble(indini(j)) + 1._dp/(2._dp*indini(j))
        	endif   
	  	enddo
	  case default
	    write(*,*)'init has an illegal value'
	    stop 

  end select

else !optimized run

 	select case (init)
	 case (0) !use input values from pfile
		  !p(1:nov)=pf(indv(1:nov)) -- pf has changed after unoptimized run
		  p(1:nov)=opf(indv(1:nov))
	case (1) !start from the result of the unoptimized run
		  p(1:nov)=pf(indv(1:nov))
	case default
	    write(*,*)'init has an illegal value'
	    stop 
	end select
endif

!write(*,'(1x,a26,1000(1x,f7.4))')'pinin        -> p(1:nov) = ',p(1:nov)

end subroutine pinin

