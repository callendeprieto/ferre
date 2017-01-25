SUBROUTINE match(a,b,n,equ)


use share, only: dp, equalize

implicit none

!in/out
integer, intent(in)					:: n		!size of input arrays a and b
real(dp),dimension(n),intent(in)	:: a,b		!input arrays
real(dp),intent(out)				:: equ      !ratio of mean/median value of arrays

!locals
real(dp)							:: mo,mf	!temporary holders of median values

equ=1.0
if (equalize == 1) then 
		equ=sum(a(1:n))/sum(b(1:n))
elseif (equalize == 2) then 
	call median(a(1:n),n,mo)
	call median(b(1:n),n,mf)
	equ=mo/mf
endif

END SUBROUTINE match

