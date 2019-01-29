
subroutine	quit(status)

! Exits graciously when status >=0
!

use, intrinsic :: iso_fortran_env, only: error_unit
use share, only: dp, snr, sffile, nfilter, lsf
use timer

implicit none
integer	:: status
real(dp):: etime

select case (status)
	case (0)
		write(error_unit,*)'quit: Run ended'
	case (1)
		write(error_unit,*)'quit: End of frd file'
	case (2)
		write(error_unit,*)'quit: End of ipf file'
	case (3)
		write(error_unit,*)'quit: End of err file'
	case (4)
		write(error_unit,*)'quit: End of wav file'
	case (10)
		write(error_unit,*)'quit: j > only_object(2)'
	case (11)
		write(error_unit,*)'quit: Error reading frd file'
	case (12)
		write(error_unit,*)'quit: Error reading ipf file'
	case (13)
		write(error_unit,*)'quit: Error reading err file'
	case (14)
		write(error_unit,*)'quit: Error reading wav file'
        case (15)
                write(error_unit,*)'quit: Error allocating memory'
	case default
		write(error_unit,*)'quit: Unrecognized status ',status
end select

call ellapsed_time(etime)

stop 1

end subroutine quit

