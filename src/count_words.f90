
 subroutine count_words(cadena,n)
!
!       Counts the number of words separated by blanks in a string
!

implicit none

character(len=*), intent(in)  :: cadena
character(len=1)				    :: a,b
integer								:: i,n	
n=0
a=cadena(1:1)
if (a.ne.' ') n=1
do i=2,len_trim(cadena)
    b=cadena(i:i)
    if(b.ne.' '.and.a.eq.' ') n=n+1 
    a=b
enddo

end subroutine count_words

