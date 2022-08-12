
subroutine qua (p,flux)


!multidimensional quadratic interpolation

use share, only: dp,ndim,indi,npix,n_p,ntot,ntimes,   &
                    f,ee,imap,                        &
		    f_format,f_access,fmtformat
				    
implicit none

!input/output 
real(dp),intent(in)	:: p(ndim)	!ndim vector of wanted pars [0-1]
real(dp),intent(out)	:: flux(npix)	!n vector of output flux
		
!local variables
integer			:: i,j,k,l	!dummy indices
real(dp)		:: t(ndim)	!input p's run 0-1, t's run 1-n_p(i)
real(dp)	 	:: wrk(npix,3**ndim) !(3**ndim) working matrix
real(dp)		:: delta	!t(indi(i))-int(t(indi(i)))
real(dp)                :: omdelta      !1-delta
real(dp)		:: delta3(3)    !array of 3 quanties derived from delta to be
					!passed to qbezier
real(dp)		:: yp		!temporary variable to hold qbezier call output
integer			:: offset(ndim) !offset=0 [-1,0,1], offset=1 [0,1,2]
integer			:: findex(3**ndim)!indices of rows in f needed for 
					  !the interpolation

!write(*,*)'in qua'

	
!	get t's from p's. The p's run between 0-1, while t's run 
!	from the min to the max value of the indices of f (1<ti<n_pi) 
t(1:ndim)=p(1:ndim)*(n_p(1:ndim)-1)+1.0_dp


offset(1:ndim)=0
if (t(1) < 2) offset(1)=1
do i=2,ndim
	if (t(i) < 2) offset(i)=1
enddo


do i=1,3**ndim	
   findex(i)=dot_product(ntimes,int(t(1:ndim))+ee(1:ndim,i)-2+offset(1:ndim))+1		
enddo


!load wrk
if (f_access == 0) then
  do i=1,3**ndim	
	wrk(1:npix,i)=f(1:npix,findex(imap(i)))		
  enddo
else
 !$omp critical
 if (f_format == 0) then 
  do i=1,3**ndim
  	read(10,fmtformat,rec=findex(imap(i))) wrk(1:npix,i)
  enddo
 else
  do i=1,3**ndim
  	read(10,rec=findex(imap(i))) wrk(1:npix,i)
  enddo
 endif
 !$omp end critical
endif


!interpolate
do i=1,ndim
	delta=t(indi(i))-int(t(indi(i)))
	omdelta=(1._dp-delta)
	delta3(1)=omdelta**2
	delta3(2)=delta**2
	delta3(3)=2._dp*delta*omdelta	
	do j=1,3**(ndim-i)
		call vqbezier(wrk(1:npix,3*j-2),wrk(1:npix,3*j-1),wrk(1:npix,3*j),& 
			delta3,wrk(1:npix,j),offset(ndim-i+1))

		!do l=1,npix
		!	call qbezier(wrk(l,3*j-2),wrk(l,3*j-1),wrk(l,3*j),& 
		!		delta3,yp,offset(ndim-i+1))
		!	wrk(l,j)=yp
		!enddo
	enddo
enddo
flux(1:npix)=wrk(1:npix,1)
			

!write(*,*)'exiting qua'

		
end subroutine qua
