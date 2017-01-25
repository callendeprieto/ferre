program kk

implicit none

integer :: i,j,icode
integer :: data(10)
character(len=10000) :: line

open(1,file='kkk')
do i=1,2
    data(:)=0.0
	read(1,*,iostat=icode) (data(j),j=1,10)
	write(*,*)i,data(1:3),icode
enddo
close(1)

end program kk