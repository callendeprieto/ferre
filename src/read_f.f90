
subroutine read_f
	
!Loading f in memory
!Checks multiple synth modules for consistency and 
!keeps track of the location of each module's fluxes,
!their resolution, type of wavelength scale (linear/log10),
!and which pixels contain photometry

!this routine is to be run serial -- no openmp critical protection provided

use share, only: dp,maxndim,maxsynth,flen,    &
		 siobuffer,liobuffer,xliobuffer,    &
		 rango, minusculo,                  &
		 constant,pcaproject,pcachi,        &
		 npca,meanspca,vpca,wpca,ff,        &
		 nelnpca,totalnpca,nvar,            &
		 n_p,ntot,ntimes,                   &
		 npix,llimits,steps, 	            & 
		 nphotpix,photpixels,		    &
		 nsynth,hs,			    &
		 ndim,nov,synthfile,fixfile,        &
		 filterfile,                        &
		 f_format,f_access,nov,cont,        &
		 f,scalef,scaled,badflux,fmtformat, &
		 winter,nfilter,lsf,npca,pcaproject,&
		 transposed,file_data19,file_data20

		 
implicit none

!locals	
real(dp), allocatable	:: 	f2(:,:)		!temporary storage for synth grid
real(dp), allocatable   ::  fixratio(:) !flux ratio correction
real(dp),allocatable    ::  record(:)	!dummy array for inquire about the record size
integer 		::  j,i,k,n_of_dim,n_p1(maxndim),npix1,npix2,status
integer			::  ii,jj,offset
integer			::  recordlength ! for unformatted synth files
integer			::  logw=0	!wavelength scale equidistant in log10?
integer			::  vacuum=0 ! wavelength scale in vacuum or std. air
integer			::  modo=-1000 ! synspec imode used
integer, allocatable	::  photpixels2(:)	!temp list of pixels with photometry
integer				::  multi = 0
integer				::  npcasynth = 0 !tracks the number of npca synth modules
character(len=flen)     ::      synthfile_internal = 'Unknown'
character(len=flen)     ::      id = 'Unknown'
character(len=flen)     ::      date = 'Unknown' 
character(len=flen)     ::      synthfile_binary
character(len=30)	::	synspec	! synspec version used
character(len=45)	::	label(maxndim),label1(maxndim)
real(dp)			:: 	llimits1(maxndim),steps1(maxndim) !phys. pars synth1
character(len=80)	::	comments1,comments2,comments3,comments4
character(len=80)	::	comments5,comments6,comments7,comments8,comments9,comments10
character(len=80)       ::      comments11,comments12,comments13,comments14 
character(len=80)       ::      speclib_vers
real(dp)			:: 	wave(2)= (/0,0/)
real(dp)       	 	::  resolution,original_sampling,continuum(4),precontinuum(4)
real(dp)			:: 	invalid_code=0.0_dp ! signals invalid entries
real(dp)			::  constant_pca=0.0_dp !tracks the value of constant for npca grids
real(dp)			::	rangef		    !range of values in f
real(dp)			::	maximo		    !max value of f
real(dp)			::	minimo		    !min value of f


namelist / synth / multi,synthfile_internal,id,date,n_of_dim,npca,n_p,npix
namelist / synth / label,llimits,steps,wave,logw,vacuum,resolution,original_sampling
namelist / synth / synspec,modo,invalid_code,constant,continuum,precontinuum
namelist / synth / transposed,file_data19,file_data20
namelist / synth / comments1,comments2,comments3,comments4,comments5,comments6
namelist / synth / comments7,comments8,comments9,comments10
namelist / synth / comments11,comments12,comments13,comments14
namelist / synth / speclib_vers


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
	if (npca(1) > 0) then
		npcasynth=1
		constant_pca=constant
		write(*,'(a39,i1,1a)')'read_f       -> NPCA grid (pcaproject=',pcaproject,')'
		if (winter == 2 .and. pcaproject == 1) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'This is an NPCA grid, but WINTER=2 (interp. library)'
			write(*,*) 'so the observed spectra should be uncompressed     '
			write(*,*) 'and pcaproject set to 0!'
			stop
		endif
		if (winter == 2 .and. pcaproject == 0) then
			if (pcachi == 0) then
				write(*,*) 'read_f: WARNING'
				write(*,*) 'This is an NPCA grid with WINTER=2 (interp. library)'
				write(*,*) 'and pcaproject=0, so the chi**2 eval. will be done on'
				write(*,*) 'uncompressed spectra'
			else
				write(*,*) 'read_f: ERROR'
				write(*,*) 'pcachi /= 0'
				write(*,*) 'This is an NPCA grid with WINTER=2 (interp. library)'
				write(*,*) 'and pcaproject=0, so the chi**2 eval. should be done on'
				write(*,*) 'uncompressed spectra and pcachi must be 0'			
			endif
		endif		
		if (winter < 2 .and. nfilter > 1) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'cannot use nfilter > 1 on an npca grid'
			write(*,*) 'unless interpolating the library (winter=2)'
			stop
		endif			
		if (lsf > 0 .and. pcachi == 1) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'pcachi = ',pcachi,' must be set to 0 when lsf > 0'
			stop
		endif
		if (filterfile .gt. ' ' .and. pcachi == 1) then
			write(*,*) 'read_f: WARNING'
			write(*,*) 'pcachi = ',pcachi,' must be set to 0 when the array'
			write(*,*) 'in filterfile is for fluxes and not PCA components'
		endif
		if(pcachi == 0. .and. pcaproject == 1) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'pcaproject = ',pcaproject,' must be set to 0 when pcachi = 0'
			stop
		endif
		if(winter == 2 .and. pcachi == 1) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'pcachi = ',pcachi,' must be set to 0 when winter = 2 with PCA grids'
			stop
		endif 	   
	else
		if (pcachi == 1) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'pcachi = ',pcachi,' must be set to 0 for non PCA grids'
			stop
		endif
		if (pcaproject == 1) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'pcaproject = ',pcaproject,' must be set to 0 for non PCA grids'
			stop
		endif		
	endif

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
		hs(1)%lambda0=wave(1)
		hs(1)%lambda1=wave(2)
		hs(1)%lws=logw
		hs(1)%lambdamin=wave(1)
		hs(1)%lambdamax=wave(1)+npix*wave(2)
                write(*,*)'lambdamin=',wave(1)
                write(*,*)'npix=',npix
                write(*,*)'wave(2)=',wave(2)
		select case (logw)
			case (1)
				hs(1)%lambdamin=10._dp**hs(1)%lambdamin
				hs(1)%lambdamax=10._dp**hs(1)%lambdamax
			case (2)
				hs(1)%lambdamin=exp(hs(1)%lambdamin)
				hs(1)%lambdamax=exp(hs(1)%lambdamax)
		end select
	endif
	
endif
do ii=1,multi
	
	read(1,nml=synth)
	
	!check that a npca constant value is not overridden
	if (npca(1) > 0 .and. abs(constant-constant_pca) > 1.e3*minusculo) then
			write(*,*) 'read_f: ERROR'
			write(*,*) 'the value of the CONSTANT keyword for an npca grid'
			write(*,*) 'has been modified'
			stop
	endif
			
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
			
			if (pcachi == 1) then
				write(*,*) 'read_f: ERROR'
				write(*,*) 'pcachi = ',pcachi,' must be set to 0 for non PCA grids'
				stop
			endif
			
		endif
		
		!keep track
		nsynth=1
		hs(1)%res=resolution
		hs(1)%pixbegin=1
		hs(1)%pixend=npix
		hs(1)%lambda0=wave(1)
		hs(1)%lambda1=wave(2)
		hs(1)%lws=logw		
		hs(1)%lambdamin=wave(1)
		hs(1)%lambdamax=wave(1)+npix*wave(2)
                write(*,*)'lambdamin=',wave(1)
                write(*,*)'npix=',npix
                write(*,*)'wave(2)=',wave(2)
		select case (logw)
			case (1)
				hs(1)%lambdamin=10._dp**hs(1)%lambdamin
				hs(1)%lambdamax=10._dp**hs(1)%lambdamax
			case (2)
				hs(1)%lambdamin=exp(hs(1)%lambdamin)
				hs(1)%lambdamax=exp(hs(1)%lambdamax)
		end select
		
	else

		!keep track
		nsynth=nsynth+1
		hs(nsynth)%res=resolution
		hs(nsynth)%pixbegin=hs(nsynth-1)%pixend+1
		hs(nsynth)%pixend=hs(nsynth-1)%pixend+npix
		hs(nsynth)%lambda0=wave(1)
		hs(nsynth)%lambda1=wave(2)
		hs(nsynth)%lws=logw	
		hs(nsynth)%lambdamin=wave(1)
		hs(nsynth)%lambdamax=wave(1)+npix*wave(2)
		select case (logw)
			case (1)
				hs(nsynth)%lambdamin=10._dp**hs(nsynth)%lambdamin
				hs(nsynth)%lambdamax=10._dp**hs(nsynth)%lambdamax
			case (2)
				hs(nsynth)%lambdamin=exp(hs(nsynth)%lambdamin)
				hs(nsynth)%lambdamax=exp(hs(nsynth)%lambdamax)
		end select
		
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
xliobuffer=liobuffer
if (npca(1) > 0) then !npca files contain means, v and w
	nelnpca=count(npca > 0)
	totalnpca=sum(npca(1:nelnpca))
	xliobuffer=25*totalnpca  !meanspca,vpca and wpca require longer rows
endif
open(1,file=synthfile(1),delim='apostrophe',recl=xliobuffer,action='read')
do ii=0,multi
	read(1,nml=synth)
enddo
npix=npix2
npix1=npix2


if (n_of_dim.ne.ndim) then
	write(*,*) 'ERROR in read_f'
	write(*,*) 'Inconsistency between the expected dimension of ',synthfile
	write(*,*) 'and the declared size ndim=',ndim
	stop
endif
if (n_of_dim.lt.nov) then
	write(*,*) 'ERROR in read_f'
	write(*,*) 'Inconsistency between the expected dimension of ',synthfile
	write(*,*) 'and number of free parameters nov=',nov
	stop
endif


allocate(ntimes(ndim))
ntot=n_p(ndim)
ntimes(ndim)=1
do j=2,ndim
        ntimes(ndim-j+1)=ntot
        ntot=ntot*n_p(ndim-j+1)
enddo


if (npca(1) > 0) then !npca files contain means, v and w

	nvar=npix/nelnpca
	allocate(meanspca(totalnpca))
	allocate(vpca(totalnpca))
	allocate(wpca(totalnpca,nvar))
	
	!read 
	read(1,*)meanspca
	read(1,*)vpca
	do i=1,nvar
		read(1,*)wpca(:,i)
	enddo
	
	write(*,'(1x,a45,4(1x,i5))')'read_f       -> totalnpca,npix,nelnpca,nvar =', &
			totalnpca,npix,nelnpca,nvar
	
	if (pcaproject == 1) then
			
		!prepare auxiliary array ff
		allocate(ff(nvar,totalnpca))
	
		k=1
		do ii=1,nelnpca
			offset=0
			if (ii > 1) offset=sum(npca(1:ii-1))
			write(*,*)'ii, offset=',ii,offset
			do i=1,npca(ii)
				do jj=1,nvar
					ff(jj,k)=wpca(k,jj)*vpca(jj+offset)
				enddo
				k=k+1
			enddo	
		enddo	
	endif
			
endif

!if synthfile access is through a file, we don't load it in memory
!but open the file for direct access and leave it open
!for the interpolation routines
if (f_access == 1) then

	!check that no more than one synthfile is used
	do j=2,maxsynth
		if (synthfile(j).gt.' ') &
			write(*,*)'read_f       -> WARNING: access-mode is direct, only 1 synthfile is used'
	enddo
	
	!no fscale or fixratio can be applied in direct-access mode
	write(*,*)'read_f       -> WARNING: access-mode is direct, scalef is not applied'
	if (fixfile(1).gt.' ') &
			write(*,*)'read_f       -> WARNING: access-mode is direct, fixfile will be ignored'

	
	!close the ascii file
	close(1)

	!open direct-access synth file
	if (f_format == 1) then 

		!allocate record and calculate recordlength
		allocate (record(npix))
		inquire (iolength=recordlength) record
    	deallocate(record)

		!compose name for binary synth file
		ii=len_trim(synthfile(1))
		synthfile_binary=synthfile(1)
		synthfile_binary=synthfile_binary(1:ii-3)
		synthfile_binary(ii-2:ii)='unf'
	
		open(10,file=synthfile_binary,access='direct',action='read', &
           form='unformatted',status='old',recl=recordlength)
           
    else

		!set recordlength and format for reading the formatted direct-access file
		recordlength=npix*12
		write(fmtformat,*) npix
		fmtformat = '(' // fmtformat(2:23) // 'ES12.5' // ')'
	
		!compose name for ascii synth file
		ii=len_trim(synthfile(1))
		synthfile_binary=synthfile(1)
		synthfile_binary=synthfile_binary(1:ii-3)
		synthfile_binary(ii-2:ii)='fmt'
	
		open(10,file=synthfile_binary,access='direct',action='read', &
           form='formatted',status='old',recl=recordlength)
            
    endif
    
    !return control to main
    return
    
endif


allocate (f(npix,ntot))	!allocate f


!if the synth file is binary we close the ascii version
!and open the binary one, then read
!otherwise, we just read the open ascii version
if (f_format == 1) then

	close(1)
	
	!calculate recordlength
	if (transposed == 0) then 
		inquire (iolength=recordlength) f(:,1)
    else
    	inquire (iolength=recordlength) f(1,:)
    endif
    
    !compose name for binary synth file
	ii=len_trim(synthfile(1))
	synthfile_binary=synthfile(1)
	synthfile_binary=synthfile_binary(1:ii-3)
	synthfile_binary(ii-2:ii)='unf'

	open(1,file=synthfile_binary,access='direct', &
           form='unformatted',status='old',recl=recordlength)

	if (transposed == 0) then
		do ii=1,ntot
			read(1,rec=ii) f(:,ii)
		enddo
	else
		do ii=1,npix
			read(1,rec=ii) f(ii,:)
		enddo		
	endif

else
	
	if (transposed == 0) then 
		do ii=1,ntot
			!i=(ii-1)/n_p(4)/n_p(3)/n_p(2)+1
			!j=(ii-(i-1)*n_p(4)*n_p(3)*n_p(2)-1)/n_p(4)/n_p(3)+1
			!k=(ii-(i-1)*n_p(4)*n_p(3)*n_p(2)-(j-1)*n_p(4)*n_p(3)-1)/n_p(4)+1
			!l=ii-(i-1)*n_p(4)*n_p(3)*n_p(2)-(j-1)*n_p(4)*n_p(3)-(k-1)*n_p(4)
			!write(*,*)i,j,k,l
			read(1,*) f(:,ii)
		enddo
	else
		do ii=1,npix
			!i=(ii-1)/n_p(4)/n_p(3)/n_p(2)+1
			!j=(ii-(i-1)*n_p(4)*n_p(3)*n_p(2)-1)/n_p(4)/n_p(3)+1
			!k=(ii-(i-1)*n_p(4)*n_p(3)*n_p(2)-(j-1)*n_p(4)*n_p(3)-1)/n_p(4)+1
			!l=ii-(i-1)*n_p(4)*n_p(3)*n_p(2)-(j-1)*n_p(4)*n_p(3)-(k-1)*n_p(4)
			!write(*,*)i,j,k,l
			read(1,*) f(ii,:)
		enddo	
	endif
	
endif

close(1)

!check if there is a fix file ratio array to be multiplied by the synth grid
if (fixfile(1).gt.' ') then

	!fixfile cannot be used when npca>0 
	if (npca(1) > 0) then 
		write(*,*) 'ERROR in read_f'
		write(*,*) 'fixfile cannot be used with PCA grids (filterfile can)'
		stop
	else
		allocate(fixratio(npix))
		open(1,file=fixfile(1),delim='apostrophe',recl=liobuffer)
		read(1,*) fixratio
		close(1)
		do ii=1,ntot
			f(:,ii)=f(:,ii)*fixratio
		enddo
		deallocate(fixratio)
	endif
endif


if (abs(resolution-0._dp) < 1e-6_dp) then 	!we deal with photometry
	allocate (photpixels(npix1))
	do j=1,npix1
		photpixels(j)=j
	enddo
	nphotpix=npix1
endif
		
		
do j=2,maxsynth

	if(synthfile(j).gt.' ') then

		write(*,*)synthfile(j)
		
		!reading header
		!use a small buffer size to read the header	
		open(1,file=synthfile(j),delim='apostrophe',recl=siobuffer)
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
			if (npca(1) > 0) then
				npcasynth=npcasynth+1
				write(*,*)'read_f       -> NPCA grid (pcaproject=',pcaproject,')'
				if (winter == 2 .and. pcaproject == 1) then
					write(*,*) 'read_f: ERROR'
					write(*,*) 'This is an NPCA grid, but WINTER=2 (interp. library)'
					write(*,*) 'so the observed spectra should be uncompressed     '
					write(*,*) 'and pcaproject set to 0!'
					stop
				endif
				if (winter == 2 .and. pcaproject == 0) then
					if (pcachi == 0) then
						write(*,*) 'read_f: WARNING'
						write(*,*) 'This is an NPCA grid with WINTER=2 (interp. library)'
						write(*,*) 'and pcaproject=0, so the chi**2 eval. will be done on'
						write(*,*) 'uncompressed spectra'
					else
						write(*,*) 'read_f: ERROR'
						write(*,*) 'pcachi /= 0'
						write(*,*) 'This is an NPCA grid with WINTER=2 (interp. library)'
						write(*,*) 'and pcaproject=0, so the chi**2 eval. should be done on'
						write(*,*) 'uncompressed spectra and pcachi must be 0'			
					endif				
				endif		
				if (winter < 2 .and. nfilter > 1) then
					write(*,*) 'read_f: ERROR'
					write(*,*) 'cannot use nfilter > 1 on an npca grid'
					write(*,*) 'unless interpolating the library (winter=2)'
					stop
				endif			
 	   			if (lsf > 0 .and. pcachi == 1) then
					write(*,*) 'read_f: ERROR'
					write(*,*) 'pcachi = ',pcachi,' must be set to 0 when lsf > 0'
					stop
				endif
				if (filterfile .gt. ' ' .and. pcachi == 1) then
					write(*,*) 'read_f: WARNING'
					write(*,*) 'pcachi = ',pcachi,' must be set to 0 when the array'
					write(*,*) 'in filterfile is for fluxes and not PCA components'
				endif
				if(pcachi == 0. .and. pcaproject == 1) then
					write(*,*) 'read_f: ERROR'
					write(*,*) 'pcaproject = ',pcaproject,' must be set to 0 when pcachi = 0'
					stop
				endif
				if(winter == 2 .and. pcachi == 1) then
					write(*,*) 'read_f: ERROR'
					write(*,*) 'pcachi = ',pcachi,' must be set to 0 when winter = 2 with PCA grids'
					stop
				endif 	
			else
				if (pcachi == 1) then
					write(*,*) 'read_f: ERROR'
					write(*,*) 'pcachi = ',pcachi,' must be set to 0 for non PCA grids'
					stop
				endif
				if (pcaproject == 1) then
					write(*,*) 'read_f: ERROR'
					write(*,*) 'pcaproject = ',pcaproject,' must be set to 0 for non PCA grids'
					stop
				endif
			endif
			
			!consistency checks for dimensions in this and previous synth modules
			do ii=1,maxndim
				if (n_p(ii).ne.n_p1(ii)) then
					write(*,*) 'ERROR in read_f'
					write(*,*)'Inconsistency between the sizes (N_P) of '
					write(*,*)'different synthfiles found when reading ',synthfile(j)
					stop
				endif
				if (abs(llimits(ii)-llimits1(ii)) > 1e-5) then
					write(*,*) 'ERROR in read_f'
					write(*,*)'Inconsistency between the contents (LLIMITS) of'
					write(*,*)'different synthfiles found when reading ',synthfile(j)
					stop
				endif
				if (abs(steps(ii)-steps1(ii)) > 1e-5) then
					write(*,*) 'ERROR in read_f'
					write(*,*)'Inconsistency between the contents (STEPS) of'
					write(*,*)'different synthfiles found when reading ',synthfile(j)
					stop
				endif
				if (label(ii).ne.label1(ii)) then
					write(*,*) 'ERROR in read_f'
					write(*,*)'Inconsistency between the labels (LABEL) of'
					write(*,*)'different synthfiles found when reading ',synthfile(j)
					stop
				endif
			enddo
			
			if (multi == 0 .and. npca(1) == 0) then
				!keep track
				nsynth=nsynth+1
				hs(nsynth)%res=resolution
				hs(nsynth)%pixbegin=hs(nsynth-1)%pixend+1
				hs(nsynth)%pixend=hs(nsynth-1)%pixend+npix
				hs(nsynth)%lambda0=wave(1)
				hs(nsynth)%lambda1=wave(2)
				hs(nsynth)%lws=logw	
				hs(nsynth)%lambdamin=wave(1)
				hs(nsynth)%lambdamax=wave(1)+npix*wave(2)
				select case (logw)
					case (1)
						hs(nsynth)%lambdamin=10._dp**hs(nsynth)%lambdamin
						hs(nsynth)%lambdamax=10._dp**hs(nsynth)%lambdamax
					case (2)
						hs(nsynth)%lambdamin=exp(hs(nsynth)%lambdamin)
						hs(nsynth)%lambdamax=exp(hs(nsynth)%lambdamax)
				end select

			endif
			
		endif

		
		do ii=1,multi
	
			read(1,nml=synth)		

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
			
			!keep track
			nsynth=nsynth+1
			hs(nsynth)%res=resolution
			hs(nsynth)%pixbegin=hs(nsynth-1)%pixend+1
			hs(nsynth)%pixend=hs(nsynth-1)%pixend+npix
			hs(nsynth)%lambda0=wave(1)
			hs(nsynth)%lambda1=wave(2)
			hs(nsynth)%lws=logw	
			hs(nsynth)%lambdamin=wave(1)
			hs(nsynth)%lambdamax=wave(1)+npix*wave(2)
			select case (logw)
				case (1)
					hs(nsynth)%lambdamin=10._dp**hs(nsynth)%lambdamin
					hs(nsynth)%lambdamax=10._dp**hs(nsynth)%lambdamax
				case (2)
					hs(nsynth)%lambdamin=exp(hs(nsynth)%lambdamin)
					hs(nsynth)%lambdamax=exp(hs(nsynth)%lambdamax)
			end select
			
			
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
		else
			open(1,file=synthfile(1),delim='apostrophe',recl=liobuffer,action='read')
		endif
		do ii=0,multi
			read(1,nml=synth)
		enddo
		npix=npix2
		
		!write(*,*)'at this point, npix1,npix=',npix1,npix
		!write(*,*)'nsynth=',nsynth

		if (n_of_dim.ne.ndim) then
			write(*,*) 'ERROR in read_f'
			write(*,*) 'Inconsistency between the expected dimension of ',synthfile
			write(*,*) 'and the declared size ndim=',ndim
			stop
		endif
		if (n_of_dim.lt.nov) then
			write(*,*) 'ERROR in read_f'
			write(*,*) 'Inconsistency between the expected dimension of ',synthfile
			write(*,*) 'and number of free parameters nov=',nov
			stop
		endif		
		
		if (npca(1) > 0) then !npca files contain means, v and w
		
			if (npcasynth > 1) then 
				write(*,*) 'ERROR in read_f'
				write(*,*) 'We are not ready yet for mixing more than one npca file!!'
				stop
			endif

			nvar=npix/nelnpca
			allocate(meanspca(totalnpca))
			allocate(vpca(totalnpca))
			allocate(wpca(totalnpca,nvar))
	
			!read 
			read(1,*)meanspca
			read(1,*)vpca
			do i=1,nvar
				read(1,*)wpca(:,i)
			enddo
	
	   write(*,'(1x,a45,4(1x,i5))')'read_f       -> totalnpca,npix,nelnpca,nvar =', &
		totalnpca,npix,nelnpca,nvar	
		
			if (pcaproject == 1) then

				!build auxiliary array ff
				allocate(ff(nvar,totalnpca))
	
				k=1
				do ii=1,totalnpca
					offset=0
					if (ii > 0) offset=sum(npca(1:ii))
					do i=1,npca(ii)
						do jj=1,nvar
							ff(jj,k)=wpca(k,jj)*vpca(jj+offset)
						enddo
						k=k+1
					enddo
				enddo	
				
			endif
			
		endif
		
		
		allocate(f2(npix1,ntot))!allocate f2
		f2=f

		
		if (abs(resolution-0._dp) < 1e-6_dp) then  	!we deal with photometry
			if (nphotpix.gt.0) then
				allocate (photpixels2(nphotpix))
				photpixels2=photpixels
				deallocate(photpixels)
				allocate (photpixels(npix+nphotpix))
				photpixels(1:nphotpix)=photpixels2
				deallocate (photpixels2,stat=status)
				do ii=1,npix
					photpixels(nphotpix+ii)=npix1+ii
				enddo
			else
				allocate (photpixels(npix))
				do ii=1,npix
					photpixels(ii)=npix1+ii
				enddo
			endif
			nphotpix=nphotpix+npix
		endif		

		npix=npix1+npix

		deallocate(f)
		allocate(f(npix,ntot))!allocate newf

		!if the synth file is binary we close the ascii version
		!and open the binary one, then read
		!otherwise, we just read the open ascii version
		if (f_format == 1) then

			close(1)
	    
	    	        !calculate recordlength
	    	        if (transposed == 0) then
				inquire (iolength=recordlength) f(npix1+1:npix,1)
			else
				inquire (iolength=recordlength) f(1,:)
			endif
			
		        !compose name for binary synth file
			ii=len_trim(synthfile(j))
			synthfile_binary=synthfile(j)
			synthfile_binary=synthfile_binary(1:ii-3)
			synthfile_binary(ii-2:ii)='unf'

			open(1,file=synthfile_binary,access='direct', &
           	   form='unformatted',status='old',recl=recordlength)

			if (transposed == 0) then
				do ii=1,ntot
					f(1:npix1,ii)=f2(:,ii)
					read(1,rec=ii) f(npix1+1:npix,ii)
				enddo
			else
				do ii=1,npix1
					f(ii,1:ntot)=f2(ii,1:ntot)
				enddo	
				do ii=npix+1,npix
					read(1,rec=ii) f(ii,1:ntot)
				enddo
			endif
			

		else

			close(1)
	    			
			open(1,file=synthfile(j),delim='apostrophe',recl=xliobuffer,action='read')
			npix2=npix
			do ii=0,multi
				read(1,nml=synth)
			enddo
			npix=npix2

			if (transposed == 0) then
				do ii=1,ntot
					f(1:npix1,ii)=f2(:,ii)
					read(1,*) f(npix1+1:npix,ii)
				enddo
			else
				do ii=1,npix1
					f(ii,1:ntot)=f2(ii,1:ntot)
				enddo	
				do ii=npix+1,npix
					read(1,*) f(ii,1:ntot)
				enddo
			endif
			
		endif
						
		close(1)
		deallocate(f2,stat=status)
		
		!check if there is a fix file ratio array to be multiplied by the synth grid
		if (fixfile(j).gt.' ') then
			
			!fixfile cannot be used when npca>0 
			if (npca(1) > 0) then 
				write(*,*) 'ERROR in read_f'
				write(*,*) 'fixfile cannot be used with PCA grids (filterfile can)'
				stop
			else
			
				allocate(fixratio(npix-npix1))
				open(1,file=fixfile(j),delim='apostrophe',recl=liobuffer)
				read(1,*)fixratio
				close(1)
				do ii=1,ntot
					f(npix1+1:npix,ii)=f(npix1+1,ii)*fixratio
				enddo
				deallocate(fixratio)
			
			endif
			
		endif
		
		npix1=npix

	endif
enddo

!write(*,*)'nsynth=',nsynth
write(*,*)''
write(*,*)'final value: npix=',npix
write(*,*)''
write(*,*)'synthfile summary'
write(*,*)'*******************************************************************'
write(*,*)'   #      R     start-pix  end-pix    start-lambda      end-lambda'
do i=1,nsynth
	write(*,'(i5,x,f9.1,i8,x,i8,x,f16.3,x,f16.3)')i,hs(i)%res, & 
          hs(i)%pixbegin,hs(i)%pixend, &
	  !hs(i)%lambda0,hs(i)%lambda1,hs(i)%lws,  &
	  hs(i)%lambdamin,hs(i)%lambdamax
enddo
write(*,*)'*******************************************************************'

!reset xliobuffer for ffile/sffile
xliobuffer=25*npix !npix data x 25 characters/datum
if (npca(1) > 0) then !npca files contain means, v and w
	nelnpca=count(npca > 0)
	totalnpca=sum(npca(1:nelnpca))
	xliobuffer=25*totalnpca  !meanspca,vpca and wpca require longer rows
endif

!scaling to avoid extreme numbers
scaled=1
if (f_access == 1 .or. npca(1) > 0 .or. nov == 0 .or. cont > 0) scaled=0
if (scaled == 1) then 
	scalef=sum(f)/real(ntot)/real(npix)
	write(*,*)'scaling ...'
	f=f/scalef
else 
	scalef=1.0_dp
endif

!define badflux value
maximo=maxval(f)
minimo=minval(f)
rangef=maximo-minimo

!basic checks for over/underflow
rango=range(f)
minusculo=tiny(f)
if (maximo > 10._dp**rango) then
	write(*,*) 'ERROR in read_f'
	write(*,*) 'over-flow: f exceeds its allowed range'
	stop
endif
if (abs(minimo) < minusculo) then
        write(*,*) 'WARNING in read_f'
        write(*,*) 'under-flow risk: f contains values < tiny'
endif

badflux=minimo-2._dp*rangef
if (badflux > -1000._dp) badflux=-1000._dp
write(*,*)'read_f       -> scalef  =',scalef
write(*,*)'read_f       -> rangef  =',rangef
write(*,*)'read_f       -> badflux =',badflux

end subroutine read_f
