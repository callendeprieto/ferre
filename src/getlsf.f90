subroutine getlsf(lsfcof,lsfarr)

!reconstucts an analytically-given lsf from its coefficients

use share, only: dp, lsf, mlsf, nlsf, dlsf

implicit none

!input/output
real(dp), intent(in)  :: lsfcof(mlsf,dlsf)
real(dp), intent(out) :: lsfarr(mlsf,nlsf)

!locals
integer				  :: i, j, offset 
real(dp)              :: x(nlsf)   ! actual array of pixels

!choose sampling
select case (lsf)
	case (3,4,13,14)
	    offset=(nlsf-1)/2
		do i=1,nlsf
			x(i)=i-1-offset
		enddo	
		!write(*,*)'nlsf=',nlsf
		!write(*,*)'mlsf=',mlsf
		!write(*,*)'x=',x
		do i=1,mlsf
			call lsf_gauss(x,nlsf,lsfcof(i,1),lsfarr(i,:))
		enddo		
	case default
	write(*,*)'getlsf: ERROR'
	write(*,*)'-- this routine only to be called for lsf=3,4,13,14'
	write(*,*)'-- not and lsf=',lsf
	stop
end select

end subroutine getlsf
