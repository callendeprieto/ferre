subroutine fsort

!sorting opfile, offile and sffile when nthreads>1
!temporary files (ending in _sorted) will be created and
!finally copied over the original ones
!

use iso_fortran_env, only: iostat_end

use share, only: dp,flen,pfile,opfile,offile,sffile, &
		siobuffer,xliobuffer,nfilter,cont

implicit none

!locals
character(len=flen)		:: opfile2,offile2,sffile2
character(len=flen) 		:: id,id2,ext
character(len=siobuffer)	:: opline 
character(len=xliobuffer)	:: ofline,sfline
integer				:: stat


ext='_sorted'
open(1,file=pfile,status='old',recl=siobuffer,action='read')
open(2,file=opfile,status='old',recl=siobuffer,action='read')
open(3,file=offile,status='old',recl=xliobuffer,action='read')
opfile2= trim(opfile) // ext
offile2= trim(offile) // ext
if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') then
	sffile2=sffile(1:len_trim(sffile)) // ext
	open(4,file=sffile,status='old',recl=xliobuffer,action='read')
	open(5,file=sffile2,recl=xliobuffer,action='write')
endif
open(7,file=opfile2,recl=siobuffer,action='write')	! output pars
open(9,file=offile2,recl=xliobuffer,action='write')	! model flux

do 
	read(1,*,iostat=stat) id
	!if (is_iostat_end(stat)) exit ! fortran2003
	if (iostat_end == stat) exit
	rewind(2)
	rewind(3)
	if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') rewind(4)
	do
		read(2,'(a)') opline
		read(3,'(a)') ofline
		if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') read(4,'(a)') sfline
		read(opline,*) id2
		if (id == id2) then
			write(7,'(a)') trim(opline)
			write(9,'(a)') trim(ofline)
		  	if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') write(5,'(a)') trim(sfline)
			exit
		endif
	enddo
enddo

close(1)
close(2)
close(3)
close(7)
close(9)
if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') then
	close(4)
	close(5)
endif


!copy sorted files to the original ones
open(2,file=opfile2,status='old',recl=siobuffer,action='read')
open(3,file=offile2,status='old',recl=xliobuffer,action='read')
if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') then
        open(4,file=sffile2,status='old',recl=xliobuffer,action='read')
        open(5,file=sffile,recl=xliobuffer,action='write')
endif
open(7,file=opfile,recl=siobuffer,action='write')       ! output pars
open(9,file=offile,recl=xliobuffer,action='write')      ! model flux

do 

  read(2,'(a)',iostat=stat) opline
  !if (is_iostat_end(stat)) exit !fortran2003
  if (iostat_end == stat) exit
  write(7,'(a)') trim(opline)

  read(3,'(a)') ofline
  write(9,'(a)') trim(ofline)

  if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') then 
    read(4,'(a)') sfline
    write(5,'(a)') trim(sfline)
  endif
  
enddo

close(2,status='delete')
close(3,status='delete')
close(7)
close(9)
if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') then
        close(4,status='delete')
        close(5)
endif


end subroutine fsort

