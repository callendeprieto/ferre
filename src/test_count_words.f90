program test_count_words

implicit none

integer			  :: n
character(len=1000) :: a

a='1 2 3 4  235.2  1211   hola'
write(*,*)'a=',a
call count_words(a,n)

write(*,*)n

end program test_count_words
