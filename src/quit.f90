
subroutine	quit(status)

! Exits graciously when status >=0
!

use share, only: dp, snr, sffile, nfilter, lsf
use timer

implicit none
integer	:: status
real(dp):: etime

select case (status)
	case (0)
		write(*,*)'quit: Run ended'
	case (1)
		write(*,*)'quit: End of frd file'
	case (2)
		write(*,*)'quit: End of ipf file'
	case (3)
		write(*,*)'quit: End of err file'
	case (4)
		write(*,*)'quit: End of wav file'
	case (10)
		write(*,*)'quit: j > only_object(2)'
	case (11)
		write(*,*)'quit: Error reading frd file'
	case (12)
		write(*,*)'quit: Error reading ipf file'
	case (13)
		write(*,*)'quit: Error reading err file'
	case (14)
		write(*,*)'quit: Error reading wav file'
	case default
		write(*,*)'quit: Unrecognized status ',status
end select

call ellapsed_time(etime)

stop

end subroutine quit

