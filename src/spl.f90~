
subroutine spl (p,flux)


!multidimensional cubic splines interpolation

use share, only: dp,ndim,indi,npix,n_p,ntot,ntimes, &
                    f,ee,imap,                      &
		    f_format,f_access,fmtformat

use booklib, only: spline_fit,spline_int
				    
implicit none

!input/output 
real(dp),intent(in)	:: p(ndim)	!ndim vector of wanted pars [0-1]
real(dp),intent(out)	:: flux(npix)	!n vector of output flux
		
!local variables
integer			:: i,j,k,l	!dummy indices
real(dp)		:: t(ndim)	!input p's run 0-1, t's run 1-n_p(i)
real(dp)	 	:: wrk(npix,4**ndim) !(4**ndim) working matrix
real(dp)            	:: delta        !t(indi(i))-int(t(indi(i)))
real(dp)		:: yp		!temporary variable to hold cbezier call output
real(dp),parameter 	:: xx0(4) = (/-1._dp,0._dp,1._dp,2._dp/)	!usual array of local abscissae
real(dp)		:: xx(4)	!local abscissae
real(dp)		:: yy(4)    !local data
real(dp)        	:: ypp(4)   !2nd derivatives
integer			:: err		!error flag for the 1D spline fitting/interp. routines
integer			:: offset(ndim) !offset=0  [-1,0,1,2] (cubic)
					  !offset=1  [0,1,2,3]  (quadratic, using  [0,1,2])
					  !offset=-1 [-2,-1,0,1](quadratic, using [-1,0,1])
integer			:: findex(4**ndim)!indices of rows in f needed for 
					  !the interpolation


!write(*,*)'in spl'

	
!	get t's from p's. The p's run between 0-1, while t's run 
!	from the min to the max value of the indices of f (1<ti<n_pi) 
t(1:ndim)=p(1:ndim)*(n_p(1:ndim)-1)+1.0_dp


offset(1:ndim)=0
do i=1,ndim
	if (t(i) < 2) offset(i)=1
	if (t(i) >= n_p(i)-1) offset(i)=-1
enddo	


do i=1,4**ndim	
   findex(i)=dot_product(ntimes,int(t(1:ndim))+ee(1:ndim,i)-2+offset(1:ndim))+1
enddo

!load wrk
if (f_access == 0) then
  do i=1,4**ndim	
	wrk(1:npix,i)=f(1:npix,findex(imap(i)))		
  enddo
else
  !$omp critical
  if (f_format == 0) then 
  	do i=1,4**ndim
  		read(10,fmtformat,rec=findex(imap(i))) wrk(1:npix,i)
  	enddo
  else
  	do i=1,4**ndim
  		read(10,rec=findex(imap(i))) wrk(1:npix,i)
  	enddo
  endif
  !$omp end critical
endif


!interpolate
do i=1,ndim
	delta=t(indi(i))-int(t(indi(i)))
	xx=xx0+real(offset(ndim-i+1))
	!if (offset(ndim-i+1) == 0) then 
		do j=1,4**(ndim-i)
			do l=1,npix
				
			   	!call cbezier(wrk(l,4*j-3),wrk(l,4*j-2), &
				!	wrk(l,4*j-1),wrk(l,4*j),delta,yp)
				
				yy(1:4)=wrk(l,4*j-3:4*j)
				call spline_fit(xx,yy,4,ypp,err)
				call spline_int(xx,yy,4,ypp,delta,yp,err)
				wrk(l,j)=yp
			enddo
		enddo
enddo
flux(1:npix)=wrk(1:npix,1)

!write(*,*)'exiting spl'

		
end subroutine spl
