
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lin (p,flux)

!multilinear interpolation

use share, only: dp, ndim,indi,npix,n_p,ntot,ntimes,     &
                            f,ee,imap,                   &
			    f_format,f_access,fmtformat
				    
implicit none

!input/output 
real(dp),intent(in)	:: p(ndim)	!ndim vector of wanted pars [0-1]
real(dp),intent(out)	:: flux(npix)	!n vector of output fluxp
		
!local variables
integer			:: i,j,k,l	!dummy indices
real(dp)		:: t(ndim)	!input p's run 0-1, t's run 1-n_p(i)
real(dp)	 	:: wrk(npix,2**ndim) !(2**ndim) working matrix	
real(dp)		:: delta	!t(indi(i))-int(t(indi(i)))
integer			:: findex(2**ndim)!indices of rows in f needed for 
					  !the interpolation


!write(*,*)'in lin'

	
!	get t's from p's. The p's run between 0-1, while t's run 
!	from the min to the max value of the indices of f (1<ti<n_pi) 
t(1:ndim)=p(1:ndim)*(n_p(1:ndim)-1)+1.0_dp

do i=1,2**ndim	
   findex(i)=dot_product(ntimes,int(t(1:ndim))+ee(1:ndim,i)-1)+1		
enddo

!load wrk
if (f_access == 0) then
  do i=1,2**ndim	
        write(*,*)'i,imap(i),findex(imap(i))=',i,imap(i),findex(imap(i))
	wrk(1:npix,i)=f(1:npix,findex(imap(i)))		
  enddo
else
  !$omp critical
  if (f_format == 0) then
  	do i=1,2**ndim
  		read(10,fmtformat,rec=findex(imap(i))) wrk(1:npix,i)
  	enddo
  else
  	do i=1,2**ndim
  		read(10,rec=findex(imap(i))) wrk(1:npix,i)
  	enddo
  endif
  !$omp end critical
endif


!interpolate
do i=1,ndim
	delta=t(indi(i))-int(t(indi(i)))
	do j=1,2**(ndim-i)
		!wrk(1:npix,j)=wrk(1:npix,2*j-1) + (wrk(1:npix,2*j)-wrk(1:npix,2*j-1))* & 
        	do concurrent (k=1:npix)
		  wrk(k,j)=wrk(k,2*j-1) + (wrk(k,2*j)-wrk(k,2*j-1))* & 
		  !(t(ndim-i+1)-int(t(ndim-i+1)))   !original
		  !(t(indi(i))-int(t(indi(i))))     !allow new interp. order 
		  delta				  !faster
                enddo
	enddo
enddo
flux(1:npix)=wrk(1:npix,1)
			

!write(*,*)'exiting lin'

		
end subroutine lin
