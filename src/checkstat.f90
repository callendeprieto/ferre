
subroutine	checkstat(stat,array_name)

! Checks the value of istat and quits throught the 'quit' function when is different from zero
!

use, intrinsic :: iso_fortran_env, only: error_unit
implicit none

integer	:: stat
character(len=30) :: array_name

if (stat /= 0) then
        write (unit=error_unit, fmt='(A, A)') &
            'Can not allocate array ', array_name
end if

end subroutine checkstat

