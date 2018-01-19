
real :: x(3) = (/1, 2 ,3/)

where (x>1) 
 x=0
 x=1
endwhere
write(*,*)x

end
