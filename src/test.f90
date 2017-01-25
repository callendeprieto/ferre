
	program test
	use booklib
	use timer
	
	implicit none
	
	real, dimension(6) 	::	x=(/0.,1.,2.,3.,4.,5./)
	real, dimension(6)	::	y=(/1.,0.,-5.,-20.,-51.,-104./)
	real, dimension(0:3)	::	c
	real                :: t
	integer			:: 	nvals=6,order=3,error,i
	
	call start_timer
	do i=1,100000
	call lsq_fit(x,y,nvals,order,c,error)
	enddo
	call ellapsed_time(t)
	write(*,*)c
	write(*,*)'t=',t
	
	end