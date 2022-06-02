subroutine msort

!sorting opfile, offile and sffile when nthreads>1
!data are held in memory for speed 
!use fsort if memory usage is an issue



use share, only: dp,flen,nobj,pfile,opfile,offile,sffile, &
		siobuffer,xliobuffer,nfilter,cont

implicit none

!locals
character(len=flen) 		:: id,id2
character(len=siobuffer)	:: opline(nobj)
character(len=xliobuffer)	:: ofline(nobj),sfline(nobj)
integer			:: i,j,stat

write(*,*)'pfile,opfile=',pfile,opfile
open(1,file=pfile,status='old',recl=siobuffer,action='read')
open(2,file=opfile,status='old',recl=siobuffer,action='read')
open(3,file=offile,status='old',recl=xliobuffer,action='read')
do i=1,nobj 
  read(2,'(a)') opline(i)
  read(3,'(a)') ofline(i)
enddo
close(2)
close(3)

open(2,file=opfile,recl=siobuffer,action='write')	! output pars
open(3,file=offile,recl=xliobuffer,action='write')	! model flux

if((nfilter > 1 .or. abs(cont) > 0) .and. sffile.gt.' ') then
	open(4,file=sffile,status='old',recl=xliobuffer,action='read')
        do i=1,nobj 
          read(4,'(a)') sfline(i)
        enddo
        close(4)
	open(4,file=sffile,recl=xliobuffer,action='write')
endif
do j=1,nobj
	read(1,*,iostat=stat) id
	do i = 1, nobj
		read(opline(i),*) id2
		if (id == id2) then
			write(2,'(a)') trim(opline(i))
			write(3,'(a)') trim(ofline(i))
		  	if((nfilter > 1 .or. abs(cont) > 0) .and. sffile.gt.' ') write(4,'(a)') trim(sfline(i))
			exit
		endif
	enddo
enddo

close(1)
close(2)
close(3)
if((nfilter > 1 .or. abs(cont) > 0) .and. sffile.gt.' ') then
	close(4)
endif




end subroutine msort

