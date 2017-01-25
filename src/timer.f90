module timer
!tracking wall time

implicit none
save

public  :: start_timer, ellapsed_time

integer,parameter       ::  r4=selected_real_kind(p=6)
integer,parameter       ::  r8=selected_real_kind(p=13)
real(r8)                ::  start_time
integer                 ::  day_of_the_month

interface ellapsed_time
	module procedure ellapsed_time_sgl
	module procedure ellapsed_time_dbl
end interface

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine start_timer 

implicit none

integer ::  time_array(8)

call date_and_time(values=time_array)

day_of_the_month = time_array(3)

start_time =  86400._r8 * time_array(3) +  & 
			  3600._r8 * time_array(5) + 60.00_r8 * time_array(6)   +   &
              time_array(7) + 0.001_r8 * time_array(8)
              
end subroutine start_timer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ellapsed_time_dbl (t)

implicit none

integer       ::  time_array(8)
real(kind=r8) ::  t

call date_and_time(values=time_array)

t = 86400._r8 * time_array(3) +  &
    3600._r8 * time_array(5) + 60.00_r8 * time_array(6)   +   &
              time_array(7) + 0.001_r8 * time_array(8)

if (time_array(3) >= day_of_the_month) then           
	t = t - start_time
	write (6, '(1x, 1a, 1f20.4,1x,a)') 'ellapsed time:  ',t,' s'
else 
	write (6,*)'error: a change in calendar month has confused the timer!'
endif

end subroutine ellapsed_time_dbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ellapsed_time_sgl (t)

implicit none

integer       ::  time_array(8)
real(kind=r4) ::  t

call date_and_time(values=time_array)

t = 86400._r8 * time_array(3) +  &
    3600._r4 * time_array(5) + 60.00_r4 * time_array(6)   +   &
               time_array(7) + 0.001_r4 * time_array(8)
           
if (time_array(3) >= day_of_the_month) then           
	t = t - start_time
	write (6, '(1x, 1a, 1f20.4,1x,a)') 'ellapsed time:  ',t,' s'
else 
	write (6,*)'error: a change in calendar month has confused the timer!'
endif
           
end subroutine ellapsed_time_sgl

end module timer

