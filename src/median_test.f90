

real, dimension(5)      ::      x=(/1.,2.,3.,4.,5./)
real					::      m

call median(x,5,m)
write(*,*)'m=',m

end