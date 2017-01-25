subroutine	convol(x,nx,lsfarr,y)

!
! 'Convolving' a flux array x with a 1D or 2D (wavelength-dependent) lsf
!

use share, only: dp,lsf,mlsf,nlsf,nsynth,hs
		 
implicit none

!input/output
real(dp),intent(in) 	:: x(nx)	  			!input data array
integer, intent(in)     :: nx 					!size of input data array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp),intent(out) 	:: y(nx)	  			!output data array


!locals
integer			:: i,j,k,ii,n,range
real(dp) 		:: tot
real(dp)     	:: y2(nx)	  					!tmp holder for output data

!write(*,*)'lsfarr in convol=', lsfarr(1,:)
!write(*,*)'x in convol=',x

range=(nlsf-1)/2
!write(*,*)'nlsf,range=',nlsf,range

!initialize y2
y2(1:nx)=0.0_dp

do j=1,nsynth		
	if (abs(hs(j)%res-0._dp) > 1e-6_dp) then	!skip photometry
		n=hs(j)%pixend-hs(j)%pixbegin+1
		do i=hs(j)%pixbegin,hs(j)%pixend

			ii=i-hs(j)%pixbegin+1
			
			!write(*,*)hs(j)%pixbegin,hs(j)%pixend

			k=1 							!1D
			if (lsf == 2 .or. lsf == 4 .or. lsf == 12 .or. lsf == 14) k=i !2D

			!write(*,*)'i,ii=',i,ii
			if (ii <= range) then
				tot=sum(lsfarr(k,2+range-ii:2*range+1))
				y2(i)=sum( lsfarr(k,2+range-ii:2*range+1)/tot * x(hs(j)%pixbegin:i+range) )
				!write(*,*)'ii <= range, lsf:',2+range-ii,2*range+1
				!write(*,*)'             x  :', hs(j)%pixbegin, i+range
			elseif (ii >= n-range) then
				tot=sum(lsfarr(k,1:range+n-ii+1 ))
				y2(i) = sum( lsfarr(k,1:range+n-ii+1 )/tot * x(i-range:hs(j)%pixend) )
				!write(*,*)'ii <= range, lsf:',1,range+n-ii+1
				!write(*,*)'             x  :', i-range,hs(j)%pixend
			else
				!write(*,*)'i,ii,range,i-range,i+range=',i,ii,range,i-range,i+range
				!write(*,*)'x(i-range:i+range)=',x(i-range:i+range)
				y2(i)=sum( lsfarr(k,:) * x(i-range:i+range) )
				!write(*,*)'lsfarr=', lsfarr(k,:)
				!write(*,*)'y2(i)=',y2(i)
			endif
		enddo	
	endif
enddo

!now safely copy to output array
y=y2

end subroutine convol
