subroutine pem(y,nel,n,yfit)

!	Break a spectrum in n pieces and normalize each 
!	by its own mean value.

use share, only: dp
implicit none

!input/output
integer, intent(in) :: nel !number of elements in y/yfit
integer, intent(in) :: n   !number of pieces
real(dp),intent(in) :: y(nel)
real(dp),intent(out):: yfit(nel)

!locals
integer             :: i
real(dp)            :: mm !mean value for pieces
integer             :: np !number of pixels per piece

np=nel/n

do i=0,n-1 
	mm=sum(y(i*np+1:np*(i+1)))/np
	yfit(i*np+1:np*(i+1))=mm
enddo
if (np*n < nel) then 
	mm=sum(y(np*n+1:nel))/(nel-np*n)
	yfit(np*i+1:nel)=mm
endif

end subroutine pem

