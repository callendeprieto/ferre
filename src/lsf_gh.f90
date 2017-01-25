subroutine lsf_gh(x,n,xcenter,par,npar,lsf)

!replicates APOGEE's lsf_gh.pro

use share, only: dp

implicit none

!input/output
real(dp), intent(in)    :: x(n)	   ! pixel location array
real(dp), intent(in)    :: par(npar)!param. array 
real(dp), intent(in)    :: xcenter ! position of the lsf center 
integer, intent(in)     :: n,npar
real(dp), intent(out)   :: lsf(n)  ! lsf array

!locals
integer,parameter       :: maxporder=5
integer					:: i,j
real(dp)				:: binsize = 0.0_dp    ! width of a pixel in x units  -- par(1)
real(dp)				:: xoffset = 0.0_dp, xcenter1 = 0.0_dp
integer, parameter	    :: Horder  = 5         ! highest allowed hermite order 
integer					:: AHorder             ! actual highest hermite order -- par(3)
integer					:: Porder(Horder+1)    ! polynomial order for variation with x
								               ! of each parameter -- par(4:Horder+4)
integer					:: cstart(Horder+2)
real(dp)				:: GHcoefs(Horder+1)   ! There are Porder[i]+1
								               ! coefficients for parameter i.  
								               ! The Hermite parameters start with H2
								               ! since we fix H1=1.
real(dp)				:: coefarr(Horder+1,maxporder+1)
real(dp)                :: cpar(npar-(Horder+4)+1)
real(dp)				:: inpar(Horder+4)



binsize=par(1)
xoffset=par(2)
AHorder = floor(par(3))
Porder(1:AHorder+1) = floor(par(4:AHorder+4))

!write(*,*)'binsize,xoffset,AHorder=',binsize,xoffset,AHorder
!write(*,*)'Porder=',Porder(1:AHorder+1)

!checks
if (binsize > 0.1) then
        write(*,*) 'lsf_gh: WARNING'
        write(*,*) 'binsize =',binsize,'> 0  '
        write(*,*)' -- but no binning will be used in building the LSF kernel'
endif
if (maxval(Porder(1:AHorder+1)) > maxporder) then
        write(*,*) 'lsf_gh: ERROR'
        write(*,*) 'maxval(Porder) = ',maxval(Porder),' > maxporder=',maxporder
        stop	
endif

cpar = par(AHorder+4+1:npar+1)
cstart(:)=0
do i=1,AHorder+1
	cstart(i+1)=cstart(i)+Porder(i)+1
enddo

!write(*,*)'cpar=',cpar
!write(*,*)'cstart=',cstart

coefarr(:,:)=0._dp
do i=1,AHorder+1
	coefarr(i,:)=0._dp
	coefarr(i,1:Porder(i)+1)= cpar(cstart(i)+1:cstart(i)+Porder(i)+1)
	!write(*,*)i,'coefarr(i,1:porder(i)+1)',coefarr(i,1:Porder(i)+1)
enddo


Xcenter1 = xcenter+xoffset

GHcoefs(:)=0._dp
do i=1,AHorder+1 
	!write(*,*)'i=',i
	GHcoefs(i) = 0.0_dp
	do j=1,Porder(i)+1
		GHcoefs(i)=GHcoefs(i)+coefarr(i,j)*Xcenter1**(j-1)
	enddo
	!write(*,*)'GHcoefs(',i,')=',GHcoefs(i)
enddo

inpar(:) = 0.0_dp
inpar(1) = 1.0_dp
inpar(2) = Xcenter 
inpar(3) = GHcoefs(1)
inpar(4) = 1.0_dp

!write(*,*)'xcenter,xoffset,xcenter1=',xcenter,xoffset,xcenter1

if (Horder > 0) then 
	do i=1,AHorder
		inpar(4+i)=GHcoefs(i+1)
	enddo
endif

!write(*,*)'n,AHorder+4=',n,AHorder+4
!write(*,*)'inpar=',inpar(1:AHorder+4)

call gausshermite(x,n,inpar(1:AHorder+4),AHorder+4,lsf)

end subroutine lsf_gh
