program ascii2bin
!
! reading an ascii synthfile and rewriting it 
! unformatted (or formatted) for fast access
!

use share, only: dp,maxndim,maxsynth,flen,  &
         siobuffer,liobuffer,xliobuffer,    &
		 rango, minusculo,                  &
		 npca,meanspca,vpca,wpca,           &
		 n_p,ntot,                          &
		 npix,llimits,steps, 		    	& 
		 nphotpix,photpixels,				&
		 nsynth,hs,							&
		 ndim,nov,synthfile,fixfile,    	&
		 f,scalef,badflux,fmtformat,        &
		 transposed,file_data19,file_data20
		 
		 
implicit none

!locals	
real(dp),allocatable ::  record(:)
integer                          ::  istat ! allocate status var
integer 			 ::  ii,j,i,n_of_dim,recordlength,npix1,npix2
integer				 ::  n_p1(maxndim)
integer				 ::  logw=0	!wavelength scale equidistant in log10?
integer			     ::  vacuum=0 ! wavelength scale in vacuum or std. air
integer				 ::	modo=-1000 ! synspec imode used
integer				 :: nelnpca = 0, totalnpca = 0
integer              :: multi = 0
character(len=flen)	 ::	synthfile_internal,id,date,synthfile_binary,synthfile_header
character(len=3)    :: style !fmt/unf
character(len=30)	 ::	synspec	! synspec version used
character(len=45)	 ::	label(maxndim),label1(maxndim)
real(dp)			 :: llimits1(maxndim),steps1(maxndim) !phys. pars synth1
character(len=80)	 ::	comments1,comments2,comments3,comments4
character(len=80)	 ::	comments5,comments6,comments7,comments8,comments9,comments10
character(len=80)       ::      comments11,comments12,comments13,comments14
character(len=80)       ::      speclib_vers
character(len=600)   :: line  ! dummy variable to copy header to hdr file
real(dp)		 	 :: wave(2)= (/0,0/)
real(dp)        	 :: resolution,original_sampling,continuum(4),precontinuum(4)
real(dp)			 :: invalid_code=0.0_dp ! signals invalid entries
real(dp)			 :: constant=0.0_dp	    ! a constant added to all data
real(dp)			 ::	rangef		    !range of values in f
real(dp)			 ::	maximo		    !max value of f
real(dp)			 ::	minimo		    !min value of f

namelist / synth / multi,synthfile_internal,id,date,n_of_dim,npca,n_p,npix
namelist / synth / label,llimits,steps,wave,logw,vacuum,resolution,original_sampling
namelist / synth / synspec,modo,invalid_code,constant,continuum,precontinuum
namelist / synth / transposed,file_data19,file_data20
namelist / synth / comments1,comments2,comments3,comments4,comments5,comments6
namelist / synth / comments7,comments8,comments9,comments10
namelist / synth / comments11,comments12,comments13,comments14
namelist / synth / speclib_vers

write(*,'(A)')'name of the f_ file to transform'
read(*,'(A)') synthfile(1)
write(*,'(A)')'fmt/unf?'
read(*,'(A)') style

j=len_trim(synthfile(1))
synthfile_binary=synthfile(1)
synthfile_binary=synthfile_binary(1:j-3)
synthfile_header=synthfile_binary
if (style.eq.'fmt') then
	synthfile_binary(j-2:j)='fmt' 
else 
	synthfile_binary(j-2:j)='unf'
endif
synthfile_header(j-2:j)='hdr'


!reading header
!use a small buffer size to read the header	
open(1,file=synthfile(1),delim='apostrophe',recl=siobuffer)
read(1,nml=synth)
npix2=0
if (multi == 0 .or. npca(1) > 0) then 
	npix2=npix
	write(*,*)'       '
	write(*,*)'read_f       -> date =',date
	write(*,*)'read_f       -> n_p  =',n_p(1:n_of_dim)
	write(*,*)'read_f       -> wave =',wave
	write(*,*)'read_f       -> id   =',id
	write(*,*)'read_f       -> label=',label(1:n_of_dim)
	write(*,*)'read_f       -> npix =',npix
	if (npca(1) > 0) write(*,*)'read_f       -> NPCA grid'

	!keep to check with other synth modules
	n_p1=n_p
	npix1=npix
	llimits1=llimits
	steps1=steps
	label1=label
	
	if (npca(1) == 0) then 
		!keep track
		nsynth=1
		hs(1)%res=resolution
		hs(1)%pixbegin=1
		hs(1)%pixend=npix
		hs(1)%lws=logw
	endif
	
endif
do ii=1,multi

	read(1,nml=synth)
			
	if (ii == 1) then 
	
		if (npca(1) == 0) then 
			write(*,*)'       '
			write(*,*)'read_f       -> date =',date
			write(*,*)'read_f       -> n_p  =',n_p(1:n_of_dim)
			write(*,*)'read_f       -> wave =',wave
			write(*,*)'read_f       -> id   =',id
			write(*,*)'read_f       -> label=',label(1:n_of_dim)

			!keep to check with other synth modules
			n_p1=n_p
			llimits1=llimits
			steps1=steps
			label1=label
			
		endif
		
		!keep track
		nsynth=1
		hs(1)%res=resolution
		hs(1)%pixbegin=1
		hs(1)%pixend=npix
		hs(1)%lws=logw		
		
	else

		!keep track
		nsynth=nsynth+1
		hs(nsynth)%res=resolution
		hs(nsynth)%pixbegin=hs(nsynth-1)%pixend+1
		hs(nsynth)%pixend=hs(nsynth-1)%pixend+npix
		hs(nsynth)%lws=logw	

		!consistency checks for dimensions in this and previous synth modules
		do i=1,maxndim
			if (n_p(i).ne.n_p1(i)) then
				write(*,*) 'ERROR in read_f'
				write(*,*)'Inconsistency between the sizes (N_P) of '
				write(*,*)'different synthfiles found when reading ',synthfile(j)
				stop
			endif
			if (abs(llimits(i)-llimits1(i)) > 1e-5) then
				write(*,*) 'ERROR in read_f'
				write(*,*)'Inconsistency between the contents (LLIMITS) of'
				write(*,*)'different synthfiles found when reading ',synthfile(j)
				stop
			endif
			if (abs(steps(i)-steps1(i)) > 1e-5) then
				write(*,*) 'ERROR in read_f'
				write(*,*)'Inconsistency between the contents (STEPS) of'
				write(*,*)'different synthfiles found when reading ',synthfile(j)
				stop
			endif
			if (label(i).ne.label1(i)) then
				write(*,*) 'ERROR in read_f'
				write(*,*)'Inconsistency between the labels (LABEL) of'
				write(*,*)'different synthfiles found when reading ',synthfile(j)
				stop
			endif
		enddo
		
		
	endif

	if (npca(1) == 0) then 		
		npix2=npix2+npix
		write(*,*)'read_f       -> npix =',npix2,'(',ii,')'
	endif
	

enddo
close(1)

!reset to a custom buffer size to read the data
liobuffer=25*npix2 !npix data x 25 characters/datum
if (npca(1) > 0) then !npca files contain means, v and w
	nelnpca=count(npca > 0)
	totalnpca=sum(npca(1:nelnpca))
	xliobuffer=25*totalnpca  !meanspca,vpca and wpca require longer rows
	open(1,file=synthfile(1),delim='apostrophe',recl=xliobuffer,action='read')
	open(2,file=synthfile_header,delim='apostrophe',recl=xliobuffer,action='write')
else
	open(1,file=synthfile(1),delim='apostrophe',recl=liobuffer,action='read')
	open(2,file=synthfile_header,delim='apostrophe',recl=liobuffer,action='write')
endif
do ii=0,multi
	line=''
	do while (line(1:1).ne.'/')
       	 read(1,'(A)') line
       	 write(2,'(A)') line
       	 line = adjustl(line)
	enddo	
enddo
npix=npix2

ntot=n_p(n_of_dim)
do j=2,n_of_dim
	ntot=ntot*n_p(n_of_dim-j+1)
enddo

if (npca(1) > 0) then !npca files contain means, v and w

	allocate(meanspca(totalnpca),stat=istat)
	call checkstat(istat,'meanspca')
	allocate(vpca(totalnpca),stat=istat)
	call checkstat(istat,'vpca')
	allocate(wpca(totalnpca,npix/nelnpca),stat=istat)
	call checkstat(istat,'wpca')
		
	!read 
	read(1,*)meanspca
	write(2,'(100000(1x,e20.10))')meanspca
	read(1,*)vpca
	write(2,'(100000(1x,e20.10))')vpca
	do j=1,npix/nelnpca
		read(1,*)wpca(:,j)
		write(2,'(100000(1x,e20.10))')wpca(:,j)
	enddo
	
endif
close(2)

!allocate record
if (transposed == 0) then 
	allocate (record(npix),stat=istat)
else
	allocate (record(ntot),stat=istat)
endif
call checkstat(istat,'record')

if (style.eq.'fmt') then 

	if (transposed == 0) then
		recordlength=npix*12
		write(*,fmtformat) npix
	else
		recordlength=ntot*12
		write(*,fmtformat) ntot
	endif
	fmtformat = '(' // fmtformat(2:23) // 'ES12.5' // ')'
	
	write(*,*)'format is ',fmtformat
	
	!open formatted file
	open(2,file=synthfile_binary,access='direct',action='write', &
	    form='formatted',status='replace',recl=recordlength)
	
	if (transposed == 0) then
		do ii=1,ntot
			read(1,*) record
			write(2,fmtformat,rec=ii) record
		enddo
	else
		do ii=1,npix
			read(1,*) record
			write(2,fmtformat,rec=ii) record
		enddo	
	endif
	
else
	!calculate recordlength
	inquire (iolength=recordlength) record

	!open binary file
	open(2,file=synthfile_binary,access='direct',action='write', &
	   	form='unformatted',status='replace',recl=recordlength)

	if (transposed == 0) then
		do ii=1,ntot
			read(1,*) record
			write(2,rec=ii) record
		enddo
	else
		do ii=1,npix
			read(1,*) record
			write(2,rec=ii) record
		enddo	
	endif
	
endif
close(1)
close(2)


end program ascii2bin
