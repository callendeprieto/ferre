
program ferre

use share
use timer
implicit none


real(dp), allocatable	   :: pf(:)		! vector of fixed parameters
real(dp), allocatable      :: pf0(:)            ! vector of params. read from pfile (in physical units) 
real(dp), allocatable	   :: spf(:)		!vector of uncertaint. all pars.
real(dp), allocatable	   :: pt(:,:)    	!storage for var. params. for individual nruns
real(dp), allocatable	   :: obs(:)		! vector of observations
real(dp)                   :: mobs          ! mean or median of obs array
real(dp), allocatable	   :: e_obs(:)		!errors in the observations
real(dp), allocatable      :: w(:)		! weights
real(dp), allocatable	   :: fit(:)		! vector of fit model fluxes (possibly wresampled)
						! (also used to hold temporarily the continuum fit to obs)
real(dp), allocatable	   :: sfit(:)		! vector of fit model fluxes 


!other (omp-shared) variables defined in data module 'share'

!locals
integer(longenough)        :: j				!counter
integer                    :: i,k,l			!+counters
integer                    :: istat			!allocate status var
integer                    :: tid,nthreads_env		   !omp
!$ integer                 :: omp_get_thread_num, omp_get_num_threads  !omp
integer	                   :: ii,jj,kk,ii1,ii2,offset,stlen   !temp variables
real(dp)                   :: ulimit,ulimit2    !temp variables
real(dp)                   :: cphot		!weight of photometry on solut.
integer		           :: ierr		!error flag for I/O
integer		           :: opterr		!error flag for optimum
real(dp)                   :: chiscale		!chi**2/sum(w_i(f_i-obs_i))^2)
real(dp)                   :: medsnr            !median S/N 
real(dp)                   :: lchi			!log10(chi**2/(nlambda1-nov+1))
real(dp), allocatable	   :: cov(:,:)		!covariance matrix (nov,nov)
real(dp), allocatable      :: ocov(:,:)	    !full cov. matrix (ndim,ndim)
real(dp), allocatable      :: bestpf(:)		!keeps track of best-fitting pars
real(dp), allocatable      :: bestspf(:)	!     uncertainties 
real(dp), allocatable      :: bestcov(:,:)  !and  covariances (abs(nrun)>1)
real(dp), allocatable      :: opf(:)		!copies of pf/spf for temp storage or
real(dp), allocatable      :: ospf(:)		!to change to absolute units and output
real(dp), allocatable      :: lambda_obs(:) !wavelength array for input data
real(dp), allocatable      :: obs_in(:)     !tmp holder for input data that need resampling
real(dp), allocatable      :: e_obs_in(:)   !tmp holder for input errors that need resampling
real(dp), allocatable	   :: obspca(:)		!tmp holder for uncompressed observations
real(dp), allocatable	   :: e_obspca(:)	!tmp holder for uncompressed errors
real(dp), allocatable      :: lsfcof1(:,:)  !holder for object-specific lsf coefs.
real(dp), allocatable      :: lsfarr1(:,:)  !holder for an object-specific lsf
real(dp)	           :: bestlchi		!keeps track of best-fitting lchi
real(dp)                   :: dum		!dummy var.
real(dp)                   :: etime             !ellapsed time (wall time, in seconds)
character(len=8)           :: date		!yyyymmdd
character(len=10)          :: time		!hhmmss.sss
character(len=5)           :: zone		!+/- hhmm offset from UTC
integer		           :: timearr(8)	!year,month,day,zone,hr,min,s,ms
character(len=flen)        :: fname		!spectrum id
character(len=flen)        :: trkfile    	!output file tracking params 
character(len=80)          :: iqheader(6)	!iq-style header holder
character(len=80)          :: iqh 		 	!singe-line header holder
character(len=48)          :: iqmargin		!margin preceding flux data
character(len=2000000)     :: dummyline     !string to hold the first line of wfile
integer		               :: status			!object status
character(len=128)          :: inputnames(100), inputfile
!=0 ok
!<0 issue with a particular object, skip object and go on
!-1: frd,  -2: ipf,   -3: err,   -10: j< only_object(1)
!>0 end of file or error, quit
!1: eof frd,2: eof ipf,3: eof err, 10: j> only_object(2)

character(len=128)           :: arg
integer		            :: ifile,nfiles

! Parse command line arguments
!  no arguments: use input.nml for single input file
!  one argument: use command line argument for single input file
!  two arguments with -l filename: use filename to provide line of input files
!  otherwise, error
if (command_argument_count() .eq. 0) then
  nfiles= 1
  inputnames(1) = 'input.nml  '
else if (command_argument_count() .eq. 1) then
  call get_command_argument(1,arg)
  nfiles= 1
  inputnames(1) = arg
else if (command_argument_count() .eq. 2) then
  call get_command_argument(1,arg)
  if ( arg .ne. '-l' ) then
    write(*,*) 'Must either provide:'
    write(*,*) '   no argument to use input.nml'
    write(*,*) '   one argument with input file name'
    write(*,*) '   two arguments with -l listname'
  endif
  call get_command_argument(2,arg)
  open(1,file=arg)
  nfiles =0
  DO
    READ(1, '(a128)', IOSTAT=istat) inputnames(nfiles+1)
    IF (istat < 0) EXIT
      nfiles = nfiles + 1
  END DO
else 
  stop 'unrecognized command line arguments'
endif

synthfile0=''
do ifile = 1, nfiles
    write(*,*) ifile, inputnames(ifile)

write(*,*)'-----------------------------------------------------------------'
write(*,'(20x,a10,10x,a12)')' f e r r e',ver
write(*,*)
call date_and_time(date,time,zone,timearr)	!find out date
write(*,'(1x,a12,i2,a1,i2,a1,i4,2x,a4,2x,i2,a1,i2,a1,f5.2)')'starting on ', &
		timearr(3),'-',timearr(2),'-',timearr(1), ' at ',timearr(5),':',    &
		timearr(6),':',real(timearr(7))+real(timearr(8))/1.e3
write(*,*)'-----------------------------------------------------------------'


!zero ttie0 and ttie
indtie(1:maxndim)=0
ttie0(1:maxndim)=0._dp
ttie(1:maxndim,1:maxndim)=0._dp

write(*,*) trim(inputnames(ifile)),len(trim(inputnames(ifile)))
call load_control(trim(inputnames(ifile)))				!read control file

!set nthreads: OMP_NUM_THREADS is used unless specified in control file
nthreads_env=1
!$omp parallel if (nthreads == 0)
!$ nthreads_env = omp_get_num_threads()
!$omp end parallel
if (nthreads == 0) nthreads=nthreads_env

				
write(*, '(a,a)') 'synthfile: ', trim(synthfile(1)), trim(synthfile0)
if (trim(synthfile(1)) .ne. trim(synthfile0)) then
call start_timer
npca(:)=0					!reset all elements of npca to 0
call read_f					!load synth grid (find out ndim,npix)
call ellapsed_time(etime)
write(*,*)'Done reading!'
synthfile0 = synthfile(1)
endif

!check for consistency between inter and n_p
!inter<min(n_p)
if (inter >= minval(n_p(1:ndim))) then
	write(*,*) 'ferre: ERROR'
	write(*,*) 'inter = ',inter,' is >= minval(n_p(1:ndim))=',minval(n_p(1:ndim))
	write(*,*) 'you need to reduce the order of the interpolations'
	stop
endif



allocate (probe(ndim),stat=istat)	!allocate variables in share module (shared for omp)
call checkstat(istat,'probe')
probe(:)=  0.2_dp			!probe > 0. ( and normally<0.5)

!setup wavelength scale, if needed
if (winter > 0) then
	if (npca(1) > 0) then 
		allocate(lambda_syn(totalnpca),stat=istat)
		call checkstat(istat,'lambda_syn')
		call setuplambda(lambda_syn,totalnpca)
	else
		allocate(lambda_syn(npix),stat=istat)
		call checkstat(istat,'lambda_syn')
		call setuplambda(lambda_syn,npix)
	endif
endif

call getee				       !load ee(ndim,2**ndim/3**ndim) 

if (nov > 0) &
  call getaa(nov,n_p(indv(1:nov))) !load aa (used in minlocus)
if (maxval(indini(1:nov)) > 0) &   
  call getuu					   !load uu (used in pinin)


!initialize random number generator (may not be needed)
call random_seed()


!need nobj?
if (nobj <= 0 .or. nthreads > 1 .or. lsf == 13 .or. lsf ==14) then ! count the number of lines in
	open(2,file=pfile,status='old',recl=siobuffer,action='read')!input file to find nobj
	write(*,*)'Counting the number of input objects...'
	j=1
	do 
		read(2,*,iostat=ierr)!fname,pf
		!write(*,*)j
		if (ierr < 0) exit
		j=j+1
	enddo
	nobj=j-1
	close(2)
endif
write(*,*)'main         -> nobj =',nobj

!need nlambda?
if (nlambda == 0 .and. winter > 0) then ! if needed, count columns in frd
	open(2,file=wfile,status='old',recl=xliobuffer,action='read')
	read(2,'(a)') dummyline
	close(2)
	!write(*,*)'dummyline(1:300)=',dummyline(1:300)
	call count_words(dummyline,nlambda)
endif
write(*,*)'main         -> nlambda =',nlambda


!setup wref and choose nlambda1 depending on wave. interp. scheme
select case (winter)
case (0)  !data and library share the same wavelength array
	nlambda1=npix
	if (npca(1) > 0 .and. pcachi == 0) nlambda1=totalnpca !unless pca library with pcachi=0
case (1)  !interpolate observations 
	nlambda1=npix
	if (lsf > 0 .and. npca(1) > 0) nlambda1=totalnpca
case (2)  !interpolate library
	nlambda1=nlambda
end select
allocate (wref(nlambda1),stat=istat)
call checkstat(istat,'wref')
wref(:) =  1.0_dp

!adopt wref from filter file if provided
if (filterfile.gt.' ') then
 	write(*,*)'main         -> filterfile = ',filterfile,' (array with ',nlambda1,' elements)'
	open(2,file=filterfile,delim='apostrophe',recl=liobuffer)
	read(2,*) wref
	close(2)
endif

!increase the weights for the photometry to count ~ half the total
!the real function should probably include a 'kind' integer match to 'dp'
!this should be checked by figuring out what 'kind' corresponds to the 'dp'
!variables, 2 is for regular doubles, but ...

if (nphotpix > 0 .and. wphot < 0._dp) then
	wphot=real(nlambda1/nphotpix) 
endif
do i=1,nphotpix
	wref(photpixels(i))=wphot
	write(*,*)'photometry:',photpixels(i),wphot
enddo

!normalize
wref=wref/sum(wref)*real(nlambda1)


!open I/O files
if (nov > 0) then !1st nov if
  open(1,file=ffile,status='old',recl=xliobuffer,action='read')		! fardo
  if (fformat == 1) then
  !skip header
	do while (iqh(1:7) .ne. '* Data ')
		read(1,'(a)',iostat=ierr) iqh
		if (abs(ierr) > 0) then
			status=1
			if (ierr > 0) status=status+10
			call quit(status)
		endif

		!write(*,*)'iqh=',iqh
		if (iqh(1:21) .eq. '* Wavelength : Start ') iqheader(1)=iqh
		if (iqh(1:21) .eq. '* Wavelength : End   ') iqheader(2)=iqh
		if (iqh(1:21) .eq. '* Wavelength : Number') iqheader(3)=iqh
		if (iqh(1:21) .eq. '* Sampling           ') iqheader(4)=iqh			
	enddo
	iqheader(5)=iqh
	
	!find out size of wavelength array
	if (nlambda == 0) read(iqh(10:24),*)nlambda
	if (winter == 0 .and. nlambda /= npix) then
		write(*,*)'ferre: ERROR'
		write(*,*)'nlambda (',nlambda,') and npix (',npix,')'
		write(*,*)'do not match, but winter = 0'
		stop
	endif
	
	!allocate waveline
	allocate (waveline(nlambda),stat=istat)
	call checkstat(istat,'waveline')
	waveline(:) = 0.0_dp		
	
	!read wavelength line
	read(1,*,iostat=ierr) iqmargin,iqh(1:1),dum,dum,dum,waveline
	if (abs(ierr) > 0) then
		status=1
		if (ierr > 0) status=status+10
		call quit(status)
	endif
  endif
else
	write(*,*)'ferre: WARNING'
	write(*,*)'Pure interpolation; no input fluxes or errors are read'
	write(*,*)'                    no output parameter file'
endif   !1st nov if

open(2,file=pfile,status='old',recl=siobuffer,action='read')	! input pars
if (nov > 0) open(3,file=opfile,recl=siobuffer,action='write')	! output pars
open(4,file=offile,recl=xliobuffer,action='write')		! model flux
if (fformat == 1) then 
	write(4,'(a)')iqheader(1:5)
	write(4,'(a9,a3,3(f5.2,7x),2000000(f12.5,1x))',advance='NO') &
				iqmargin,' : ',dum,dum,dum,waveline
	write(4,'(f12.1)') -2.
endif
if (snr == -1._dp .and. nov > 0) open(5,file=erfile,status='old',recl=liobuffer,action='read') ! err
if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') then
	open(7,file=sffile,recl=xliobuffer,action='write')  ! smoothed obs flux
	if (fformat == 1) then
		write(7,'(a9,a3,3(f5.2,7x),2000000(f12.5,1x))',advance='NO') &
				iqmargin, ' : ',dum,dum,dum,waveline
		write(7,'(f12.1)') -2.
	endif
endif
if (wfile.gt.' ' .and. winter > 0) open(8,file=wfile,status='old',recl=liobuffer,action='read')! wav

!general lsf
select case (lsf)
	case (0,1,3,11,13)
		mlsf=1
		dlsf=1  !relevant when lsf=3,13
	case (2,4,12,14)
		mlsf=npix
                if (npca(1) > 0) mlsf=totalnpca
		dlsf=1  !relevant when lsf=4,14
	case default
		write(*,*)'ferre: ERROR'
		write(*,*)'-- lsf must be 0,1,2,3,4,11,12, 13 or 14 '
		write(*,*)'-- this should have been caught in load_control!'
		stop
end select

if (lsf > 0 .and. lsf < 10) then !lsf is the same for all objects
	if (lsf < 3) then !is read from file
		allocate(lsfarr(mlsf,nlsf),stat=istat)
		call checkstat(istat,'lsfarr')
		lsfarr(:,:)=0.0_dp
		open(9,file=lsffile,status='old',recl=xliobuffer,action='read')
		read(9,*,iostat=ierr) lsfarr
		close(9)
	else  ! we read lsfcof from file
		allocate(lsfcof(mlsf,dlsf),stat=istat)
		call checkstat(istat,'lsfcof')
		lsfcof(:,:)=0.0_dp
		open(9,file=lsffile,status='old',recl=xliobuffer,action='read')
		read(9,*,iostat=ierr) lsfcof
		write(*,*)'main         -> lsfcof=',lsfcof
		close(9)
		nlsf=ceiling(1.4*maxval(lsfcof(:,1)))*2+1 !set nlsf
		write(*,*)'main         -> nlsf = ',nlsf 
		allocate(lsfarr(mlsf,nlsf),stat=istat)
		call checkstat(istat,'lsfarr')
		lsfarr(:,:)=0.0_dp		
		!and calculate lsfarr
		call getlsf(lsfcof,lsfarr)
	endif
endif

!if lsfarr is parametrized and depends on the object, we
!need to scan the entire set of coefficients to decide on nlsf
if (lsf == 13 .or. lsf == 14) then
	allocate(lsfcof(mlsf,dlsf),stat=istat)
	call checkstat(istat,'lsfcof')
	lsfcof(:,:)=0.0_dp
	open(9,file=lsffile,status='old',recl=xliobuffer,action='read')
	j=0
	do i=1,nobj
		read(9,*,iostat=ierr) lsfcof
		nlsf=ceiling(1.4*maxval(lsfcof(:,1)))*2+1 !set nlsf
		if (nlsf > j) j=nlsf
	enddo
	nlsf=j
	write(*,*)'main         -> nlsf = ',nlsf
	close(9)
endif


write(*,*)'about to enter parallel region!'


!$omp parallel num_threads(nthreads)                            	&
!$omp	default(none)                                    	        &
!$omp   private(pf,pf0,spf,pt,obs,mobs,e_obs,w,fit,sfit,                &
!$omp		 j,i,k,l,    						&
!$omp   	 istat,                               		        & 
!$omp	         tid,nthreads_env,					&
!$omp		 ii,ii1,ii2,jj,kk,offset,stlen,                       	&
!$omp            ulimit,ulimit2,cphot,ierr,opterr,	                &
!$omp		 chiscale,medsnr,lchi,	  				&
!$omp		 cov,ocov,bestpf,bestspf,opf,ospf,bestlchi,bestcov, 	&
!$omp            lambda_obs,obs_in,e_obs_in,                        	&
!$omp		 obspca,e_obspca,                                       &
!$omp            lsfcof1,lsfarr1,                                   	&
!$omp		 dum,etime,date,time,zone,timearr,			&
!$omp		 fname,trkfile,iqheader,iqh,iqmargin,status)		&
!$omp	shared(npca,meanspca,vpca,wpca,ff,                   		&
!$omp            totalnpca,nelnpca,nvar,                            	&
!$omp            pcaproject,pcachi,constant,                        	&
!$omp            lsf,mlsf,nlsf,dlsf,lsfcof,lsfarr,                 	&
!$omp            n_p,npix,llimits,steps,                    	    	&
!$omp            nphotpix,photpixels,scalef,scaled,    		    	&
!$omp	         siobuffer,liobuffer,xliobuffer,   	  	    	&
!$omp		 rango,minusculo,                   	            	&
!$omp	         nsynth,hs,					    	&
!$omp            ndim,nov,indv,indini,                                  &
!$omp            ntie,indtie,typetie,ttie0,ttie,                        &
!$omp            nlambda,nobj,                 	                        &
!$omp		 synthfile,fixfile,filterfile,pfile,ffile,erfile,   	&
!$omp		 opfile,offile,sffile,lsffile,wfile,       	    	&
!$omp            f_format,f_access,fformat,snr,only_object,         	&
!$omp            ycutoff,wphot,balance,optimize,impact,             	&
!$omp            mforce,chiout,trkout,cont,ncont,obscont,rejectcont,	&
!$omp            nfilter,init,nruns,errbar,covprint,indi,    	&
!$omp            inter,mono,algor,scope,stopcr,simp,                	&
!$omp            nlambda1,winter,                                   	&
!$omp            probe,                                             	&
!$omp            f,aa,ee,uu,imap,                                  	&
!$omp            badflux,fmtformat,                                 	&
!$omp            wref,waveline,lambda_syn,                          	&
!$omp            nthreads,                                          	&
!$omp		 start_time,day_of_the_month) ! module timer 

!confirm nthreads in actual parallel loop
!$omp master
!$ nthreads=omp_get_num_threads()
write(*,*)'main         -> nthreads =',nthreads
!$omp end master

!threads identify themselves
tid=1
!$ tid=omp_get_thread_num()+1
write(*,*)'main         -> tid =',tid

!allocate locals that are private for omp
allocate (pf(ndim),stat=istat)
call checkstat(istat,'pf')
allocate (opf(ndim),stat=istat)
call checkstat(istat,'opf')
allocate (pf0(ndim),stat=istat)
call checkstat(istat,'pf0')
allocate (spf(ndim),stat=istat)
call checkstat(istat,'spf')
if (errbar == -2 .and. nruns > 1) then
	allocate(pt(nruns,nov),stat=istat)
	call checkstat(istat,'pt')
endif
allocate (ospf(ndim),stat=istat)
call checkstat(istat,'ospf')
allocate (bestpf(ndim),stat=istat)
call checkstat(istat,'bestpf')
allocate (bestspf(ndim),stat=istat)
call checkstat(istat,'bestspf')
allocate (cov(nov,nov),stat=istat)
call checkstat(istat,'cov')
allocate (bestcov(nov,nov),stat=istat)
call checkstat(istat,'bestcov')
if (covprint == 1) then
	allocate (ocov(ndim,ndim),stat=istat)
	call checkstat(istat,'ocov')
endif

!allocate w,obs,e_obs,fit, and obs_in/lambda_obs when needed
select case (winter)
case (0)  !data and library share the same wavelength array
	!nlambda1=npix
case (1)  !interpolate observations 
	!nlambda1=totalnpca when lsf>0 and npca grid
	!nlambda1=npix otherwise
	allocate(lambda_obs(nlambda),stat=istat)
	call checkstat(istat,'lambda_obs')
	allocate(obs_in(nlambda),stat=istat)
	call checkstat(istat,'obs_in')
	allocate(e_obs_in(nlambda),stat=istat)
	call checkstat(istat,'e_obs_in')
case (2)  !interpolate library
	!nlambda1=nlambda
	allocate(lambda_obs(nlambda),stat=istat)
	call checkstat(istat,'lambda_obs')
end select
	
allocate (w(nlambda1),stat=istat)			
call checkstat(istat,'w')
allocate (obs(nlambda1),stat=istat)		
call checkstat(istat,'obs')
allocate (e_obs(nlambda1),stat=istat)		
call checkstat(istat,'e_obs')
allocate (fit(nlambda1),stat=istat)
call checkstat(istat,'fit')
allocate (sfit(npix),stat=istat)
call checkstat(istat,'sfit')

if (npca(1) > 0) then 
	allocate(obspca(totalnpca),stat=istat)
	call checkstat(istat,'obspca')
	allocate(e_obspca(totalnpca),stat=istat)
	call checkstat(istat,'e_obspca')
endif

!allocate object-specific lsfarr1 when needed
allocate(lsfarr1(mlsf,nlsf),stat=istat)
call checkstat(istat,'lsfarr1')
allocate(lsfcof1(mlsf,dlsf),stat=istat)
call checkstat(istat,'lsfcof1')
if (lsf > 10) open(9,file=lsffile,status='old',recl=xliobuffer,action='read')

!processing
status=0
!$omp do
do j=1,nobj


	write(*,*)'next object #',j
	
	tid=1
	!$ tid=omp_get_thread_num()+1

	!reading input data
	!$omp critical
	
	!reading1 -- spectra id and parameters
	pf(:)=-1.0_dp
	pf0(:)=-1.0_dp
  	spf(:)=-1.0_dp
	!read(2,'(1x,a30,40(1x,f8.5))',err=200,end=100) fname,pf,spf
	read(2,*,iostat=ierr) fname,pf
        if (trkout /= 0)  then 
		stlen=len_trim(fname)
		trkfile=''
     		trkfile(1:stlen)=fname
		trkfile(stlen+1:stlen+4)='.trk'
		open(11,file=trkfile,status='unknown',recl=siobuffer,action='write') ! track parameters
		if (trkout < 0) then
		  trkfile(stlen+1:stlen+4)='.frk'
                  open(12,file=trkfile,status='unknown',recl=xliobuffer,action='write') ! track fitting residuals
                endif 
	endif

	pf0(1:ndim)=pf(1:ndim) ! pf0 will retain original values, 
			       ! in physical units, for type 1 ties
	
	if (abs(ierr) > 0) then
	status=2
	if (ierr > 0) status=status+10
		call quit(status)
	endif
	
	!reading2 -- lsf 
	lsfarr1(:,:)=0.0_dp
	if (lsf > 0) then
		select case (lsf)
		case (1,2,3,4)  ! same lsf for all objects
			lsfarr1=lsfarr
		case (11,12)! a fresh lsf for each object from lsffile
			read(9,*)lsfarr1
			!check that lsf has been normalized!
		case (13,14) ! a fresh lsf for each object, reconstructed from lsfcof1
			read(9,*)lsfcof1
			call getlsf(lsfcof1,lsfarr1)
		case default
			write(*,*)'ferre: ERROR'
			write(*,*)'-- lsf must be 0,1,2,3,4,11,12, 13 or 14 '
			write(*,*)'-- this should have been caught in load_control!'
			stop
		end select
	endif	

	!initializing e_obs
	e_obs(:)=1._dp
	
	!write(*,*)'lsfarr1 in ferre=',lsfarr1(1,:)
	  
        if (nov > 0) then 	!2nd nov if
	
	  !reading3 -- observed spectrum and wavelengths
 	  obs(:)=0.0_dp
	
 	  select case(winter)
 	  case (0)  !no interpolation in wavelength needed
 	    if (npca(1) > 0 .and. pcaproject == 1) then 
	 	  obspca(:)=0.0_dp
		  if (fformat == 1) then
			read(1,*,iostat=ierr) iqmargin,iqh,dum,dum,dum,obspca(:)
			write(*,*) iqmargin !,obs(:)			
		  else
			read(1,*,iostat=ierr) obspca(:)	
		  endif		
            else
                  if (fformat == 1) then
			read(1,*,iostat=ierr) iqmargin,iqh,dum,dum,dum,obs(:)
			write(*,*) iqmargin !,obs(:)			
		  else
			read(1,*,iostat=ierr) obs(:)		
		  endif
	    endif
 	  case (1)   !interpolate observations
 	    if (npca(1) > 0 .and. pcaproject == 1) then
 	      obspca(:)=0.0_dp
 	      if (fformat == 1) then 
 	    	read(1,*,iostat=ierr) iqmargin,iqh,dum,dum,dum,obs_in(:)
 	    	lambda_obs(1:nlambda)=waveline(1:nlambda)
 	      else 
 	    	read(1,*,iostat=ierr) obs_in(:)
 	    	read(8,*,iostat=ierr) lambda_obs(:)
 	      endif
 	      call wresample(lambda_obs,obs_in,nlambda,lambda_syn,obspca,totalnpca)
 	    else
 	      if (fformat == 1) then 
 	    	read(1,*,iostat=ierr) iqmargin,iqh,dum,dum,dum,obs_in(:)
 	    	lambda_obs(1:nlambda)=waveline(1:nlambda)
 	      else 
 	    	read(1,*,iostat=ierr) obs_in(:)
 	    	read(8,*,iostat=ierr) lambda_obs(:)
 	      endif
 	      if (pcaproject == 0 .and. npca(1) > 0 .and. lsf > 0) then
 	      	call wresample(lambda_obs,obs_in,nlambda,lambda_syn,obs,totalnpca)
 	      else
 	      	call wresample(lambda_obs,obs_in,nlambda,lambda_syn,obs,npix)	      	
 	      endif
 	    endif
 	  case (2)   !interpolate library
 	      if (fformat == 1) then 
 	    	  read(1,*,iostat=ierr) iqmargin,iqh,dum,dum,dum,obs(:)
 	     	  lambda_obs(1:nlambda)=waveline(1:nlambda)
 	      else 
 	    	  read(1,*,iostat=ierr) obs(:)
 	    	  read(8,*,iostat=ierr) lambda_obs(:)
 	      endif	    
 	  case default
 	    write(*,*)'ferre: ERROR'
 	    write(*,*)'winter has an ilegal value of ',winter
 	    stop
	  end select
	  
	  if (abs(ierr) > 0) then
		status=1
		if (ierr > 0) status=status+10
		call quit(status)
	  endif      
	  
	  !reading4 -- spectrum errors
	  if (snr == -1._dp) then !read flux errors
	  	e_obs(:)=1.0_dp	
	  	select case (winter)
	  	case (0)
		  	if (npca(1) > 0 .and. pcaproject == 1) then
				read(5,*,iostat=ierr) e_obspca
			else 
				read(5,*,iostat=ierr) e_obs
			endif		
	  	case (1) 
		  	if (npca(1) > 0 .and. pcaproject == 1) then
				read(5,*,iostat=ierr) e_obs_in
				call wresample(lambda_obs,e_obs_in,nlambda,lambda_syn,e_obspca,totalnpca)
			else 
				read(5,*,iostat=ierr) e_obs_in
				call wresample(lambda_obs,e_obs_in,nlambda,lambda_syn,e_obs,npix)			
			endif				  	
	  	case (2)
		  	if (npca(1) > 0 .and. pcaproject == 1) then
				read(5,*,iostat=ierr) e_obs
			else 
				read(5,*,iostat=ierr) e_obs
			endif	
		end select
		if (abs(ierr) > 0) then
			status=3
			if (ierr > 0) status=status+10
			call quit(status)
	  	endif
	  endif	  
    
	  !make sure we do not have errors=0 (used on the error calcs.)
	  if (minval(e_obs) <= tiny(e_obs)) then
	        write(*,*) 'ferre: WARNING'
        	write(*,*) 'e_obs contains values < tiny'
		write(*,*) 'which are reset to ',maxval(e_obs)*1.e6_dp+1._dp
		where (e_obs < tiny(e_obs)) e_obs=maxval(e_obs)*1.e6_dp+1._dp
	  endif

	  
	else    !even when nov =0 we may want to read the observed/targeted wavelengths for interpolation
	  if (wfile.gt.' ' .and. winter > 0) then 
 	      if (fformat == 1) then 
 	    	  read(1,*,iostat=ierr) iqmargin,iqh,dum,dum,dum,obs(:)
 	     	  lambda_obs(1:nlambda)=waveline(1:nlambda)
 	      else 
 	    	  read(1,*,iostat=ierr) obs(:)
 	    	  read(8,*,iostat=ierr) lambda_obs(:)
 	      endif	   	
 	   endif
	endif	! 2nd nov if

	!$omp end critical

        !write(*,*) 'status,tid,pf=',status,tid,pf
	!check status
	if (status == 0) then  !check1 status	
	  if (nov > 0) then    !3rd nov if
		

	  	!npca projection
	  	if (npca(1) > 0 .and. pcaproject == 1) then
	  		obspca = obspca - meanspca
			do ii=1,nelnpca
		 	  offset=0
		  	 if (ii > 1) offset=sum(npca(1:ii-1))
				do jj=1,nvar 
					obs((ii-1)*nvar+jj)=0.0_dp
					e_obs((ii-1)*nvar+jj)=0.0_dp
					do i=1,npca(ii)
						kk=i+offset
						obs((ii-1)*nvar+jj) = obs((ii-1)*nvar+jj) + obspca(kk)*ff(jj,kk)
						e_obs((ii-1)*nvar+jj) = e_obs((ii-1)*nvar+jj) + &
							e_obspca(kk)**2*ff(jj,kk)**2
					enddo
					e_obs((ii-1)*nvar+jj)=sqrt(e_obs((ii-1)*nvar+jj))
				enddo
			enddo
			!write(*,*)obspca + meanspca
			obs = obs + constant  	
	  	endif
	  	

	  	!apply scaling
	  	if (scaled == 1) obs=obs/scalef

	  	
		!determine mean/median for obs when mforce>0
		mobs=0.0_dp
		if (mforce > 0) then
		    if (mforce == 1) then
	 	 	   mobs=sum(obs(1:nlambda1))/nlambda1
		    else 
		 	   call median(obs,nlambda1,mobs)
		    endif
		endif	  

			
	  	!assign reference weights
	  	w=wref
	  
	 	!balance weights if requested
	  	if (balance.eq.1) then 
	  			where (abs(obs) > tiny(obs)) w=1._dp/obs**2
	  	endif
	  

		!restrict wavelength to use in chi**2 if ycutoff>0
		where (obs < ycutoff) w = 0.0_dp
		if (pcachi == 0) where (obs <= 0.0_dp) w = 0.0_dp
		if (sum(abs(w)) >  0.0_dp) then 
			!renormalize
			w=w/sum(w)*real(nlambda1)
	  	else
			status=-1		
	  	endif	


	    if (snr == -1._dp) then !use flux errors
		
			!apply scaling
			if (scaled == 1) e_obs=e_obs/scalef
			
		
			!boxcar smoothing
			if (nfilter >  1) call smooth2(obs,e_obs,nlambda1,nfilter)			

	
	   	 	!
	    	!chi**2 = sum_i w(i)/e_obs(i)^2 (flux(i)-obs(i))^2 
	    	!	= chiscale * sum_i v(i) (flux(i)-obs(i))^2
	    	!	
	    	!	sum_i w(i)=npix 
	    	!	chiscale=sum_j w(j)/e_obs(j)^2/npix
	    	!	v(i)= w(i)/e_obs(i)^2/chiscale
	    	!	sum_i v(i)=npix
	    	!
	    	!	The original w weights take care of bad points 
	    	!	(obs < ycutoff .or. obs <= 0.0_dp)
	    	!	and balance between different datasets 
	    	!	(e.g. photom. vs. spec, as specified in wref).
	    	!	We derive v(i), which are introduced into w
	    	!	and save chiscale to calculate later the proper chi**2.
	    	!
	    	
	
			!set e_obs to 1 when w=0 to ensure those data are not considered in chi2	
!			where (w == 0.0_dp) e_obs=1._dp
			

			chiscale=sum(w/e_obs**2)/real(nlambda1)
			
			if (abs(chiscale) >  0.0_dp) then 
			    w=w/e_obs**2/chiscale
			else
			    status=-3
                            write(*,*)'This object appears to have all data equal to zero. It will be skipped.'
			endif
		
			!write(*,*)'weights=',w(1:14)
			!write(*,*)'chiscale=',chiscale
			!write(*,*)'total weights phot/spec:',sum(w(1:4)),sum(w(5:701))


	    else
	
	    	!	When all e_obs(i) are equal or unknown (and therefore snr is
	    	!	provided in the control file and > 0) then
	    	!chi**2 = sum_i w(i)/e_obs(i)^2 (flux(i)-obs(i))^2 
	   	 	!	= snr**2 sum_i w(i) (flux(i)-obs(i))^2
	    	!	and therefore chiscale=snr^2 and v(i)=w(i).	
	
			chiscale=snr**2
		
			!boxcar smoothing
			if (nfilter >  1) call smooth1(obs,nlambda1,nfilter)
		
			!load e_obs
			e_obs=abs(obs/snr)
		
	    endif

	
	    !continuum normalization
	    if (cont>0 .and. obscont /= 0) then 
	      call continuum(obs,lambda_obs,e_obs,fit,nlambda1,cont,ncont,rejectcont)
              where (fit /= 0._dp)
	          obs=obs/fit
	          e_obs=e_obs/fit
              endwhere
	    endif


	  endif  !3rd nov if
	
	
	  !from physical to normalized ([0-1]) units
	  call normal(pf)

	  !store parameter values from input file (may be required for optimized fits)
	  opf=pf
	  
	  if (nov > 0) then
	    !write(*,*)'Done reading input flux #',j
	  	!write(*,'(A,1000(1X,I3))')' indv=',indv(1:nov)
	  endif
	  
	  !checking that all the fixed input parameters and also the variable 
	  !parameters initialized from pfile (init=0) are within [0,1]
	  if (any(pf < 0.0_dp .or. pf >= 1.0_dp)) then
		do i=1,ndim
			if (pf(i) < 0.0_dp .or. pf(i) >= 1.0_dp) then 
				status=-2	
				do k=1,nov
				  if (indv(k) == i .and. init>0) status=0
				enddo
			endif
		enddo
		if (status == -2) then
			write(*,*) 'Input parameters are outside the range'
		endif
	  endif
	
	  !only_object keyword to skip objects
	  if (j < only_object(1)) status=-10

	endif			!check1 status	
	
	bestpf(:)=-1.0_dp
	bestspf(:)=-1.0_dp
	bestlchi=1.e10_dp
	bestcov(:,:)=0.0_dp

	do k=1,abs(nruns)	!loop on nruns/object
	
	if (status == 0 .and. nov > 0) then   !check2 status and 4th nov if

	  call getmin(algor,k,0,fname,chiscale, & 
		      w,pf,pf0,opf,obs,lambda_obs,e_obs,mobs,lsfarr1,& 
		      spf,lchi,cov)  

	  !keep track of results when using nrunsigma
	  if (errbar == -2 .and. nruns>1) then 
		pt(k,1:nov)=pf(indv(1:nov))
          	if (k >= nruns) then 
	  	!calculate cov. and errors from the nruns solutions
	  		call nrunsigma(pt,spf,cov)
		endif
	  endif
	  	 
      	  do i=1,nov
            if (abs(spf(indv(i))).ge.999) spf(indv(i))=-1.0_dp
      	  enddo            

	  WRITE (*,*)
	  WRITE (*,*)j,fname
	  WRITE (*,'(A4,1X,4(f8.4,1X))')  ' SOL',pf(1:ndim)
	  WRITE (*,'(A4,1X,4(f8.4,1X))')  ' ERR',spf(1:ndim)


	  !optimize weights and rerun
	  if(optimize == 1 .and. lchi < 4.0_dp) then
	  
	      call optimum(obs,lambda_obs,e_obs,mobs,lsfarr1,pf,w,opterr)
		
	      if (opterr.eq.0) then 		!check opterr
		
	    !restrict wavelength to use in chi**2 if ycutoff>0
		!(the extra-weight to phot pixels is taken care of 
		! inside optimum)
		where (obs < ycutoff) w = 0.0_dp
		if (npca(1) == 0) where (obs <= 0.0_dp) w = 0.0_dp
		if (pcachi == 0) where (obs <= 0.0_dp) w = 0.0_dp
		!renormalize
		if (sum(abs(w)) >  0.0_dp) then 
			!renormalize
			w=w/sum(w)*real(nlambda1)
	  	else
			status=-1		
	  	endif			
		
		if (snr == -1._dp) then !use flux errors
			!set e_obs to 1 when w=0 (avoiding division by 0)	
			where (w == 0.0_dp) e_obs=1._dp
			chiscale=sum(w/e_obs**2)/real(nlambda1)
			if (abs(chiscale) >  0.0_dp) then 
				w=w/e_obs**2/chiscale
			else
				status=-3
			endif			
		else
			chiscale=snr**2
		endif		

	  	call getmin(algor,k,1,fname, chiscale, & 
			    w,pf,pf0,opf,obs,lambda_obs,e_obs,mobs,lsfarr1,& 
			    spf,lchi,cov)  

	  	!keep track of results when using nrunsigma
	  	if (errbar == -2 .and. nruns>1) then 
			pt(k,1:nov)=pf(indv(1:nov))
          		if (k >= nruns) then 
	  		!calculate cov. and errors from the nruns solutions
	  			call nrunsigma(pt,spf,cov)
			endif
	  	endif
	  	 
      	  	do i=1,nov
          	  if (abs(spf(indv(i))).ge.999) spf(indv(i))=-1.0_dp
      	  	enddo            	

		WRITE (*,*)j,fname
		WRITE (*,'(A4,1X,4(f8.4,1X))')  ' SOL',pf(1:ndim)
		WRITE (*,'(A4,1X,4(f8.4,1X))')  ' ERR',spf(1:ndim)
		
	      endif	!check opterr				
	  endif ! optimization	
	 
	endif ! check2 status and 4th nov if
	
	if (lchi<bestlchi) then
		bestpf=pf
		bestspf=spf
		bestlchi=lchi
		bestcov=cov
	endif
			


	!writing output file(s)
	if (status > -10 .and. k>= nruns) then 

		if (nruns>1 .and. k==nruns) then 
			pf(1:ndim)=bestpf(1:ndim)
			lchi=bestlchi
			if (errbar /= -2) then 
				spf(1:ndim)=bestspf(1:ndim)
				cov(1:nov,1:nov)=bestcov(1:nov,1:nov)
			endif
		endif

                !getting and writing model fluxes
                call flx(pf,lambda_obs,e_obs,mobs,lsfarr1,fit)

                !$omp critical			
		!writing smoothed/normalized observed fluxes
		if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') then
			if (fformat == 0) then
			    write(7,'(2000000(es12.5,1x))') obs(1:nlambda1)*scalef
			else if (fformat == 1) then
			    write(7,'(a9,a3,3(f5.2,8x),2000000(es12.5,1x))', &
			    advance='NO') iqmargin,' : ',dum,dum,dum,obs(1:nlambda1)*scalef
				write(7,'(F12.1)') -2.
			endif
		endif


	  	!npca (de)projection if necessary
	  	if (npca(1) > 0 .and. pcaproject == 1 .and. nlambda1 < totalnpca) then
			call decompress(fit,obspca)					
			if (fformat == 0) then
			   	write(4,'(2000000(es12.5,1x))') obspca(1:totalnpca)
			else if (fformat == 1) then
			   	write(4,'(a9,a3,3(f5.2,7x),2000000(es12.5,1x))', &
			   	advance='NO')iqmargin,' : ',dum,dum,dum,obspca(1:totalnpca)
			 		write(4,'(F12.1)') -2.	
			endif
		else		
			if (fformat == 0) then
			  	write(4,'(2000000(es12.5,1x))') fit(1:nlambda1)*scalef
			else if (fformat == 1) then
			   	write(4,'(a9,a3,3(f5.2,7x),2000000(es12.5,1x))', &
			   	advance='NO')iqmargin,' : ',dum,dum,dum,fit(1:nlambda1)*scalef
			 		write(4,'(F12.1)') -2.	
			endif
		endif
		
	
		!output chi**2 surface
		if (chiout.gt.0) call chisurf(w,obs,lambda_obs,fname)

		!from  normalized to physical units
		do i=1,ndim
		  ulimit=llimits(i)+steps(i)*(n_p(i)-1)
		  opf(i)=pf(i)*(ulimit-llimits(i))+llimits(i)
		  ospf(i)=spf(i)*(ulimit-llimits(i))
		enddo	

		if (covprint == 1) then
		  ocov(:,:)=0.0_dp
		  do i=1,nov
			ii1=indv(i)
			ulimit=llimits(ii1)+steps(ii1)*(n_p(ii1)-1)
			do l=1,nov
				ii2=indv(l)
				ulimit2=llimits(ii2)+steps(ii2)*(n_p(ii2)-1)
				ocov(ii2,ii1)=cov(l,i)*(ulimit-llimits(ii1))*(ulimit2-llimits(ii2))
			enddo
		  enddo
		endif

		where (ospf < 0. .or. ospf > 99999.) 	ospf=-999.999_dp
		where (opf > 99999.)		  	 opf=-999.999_dp
		if (status /=0) then
		  ospf=-999.999_dp
		  opf=-999.999_dp
		  chiscale=1.0_dp
		  lchi=-999.999_dp 
		endif
	
		!contribution made by the photometry
		!to the final solution (weight of spectroscopy  = 1.-cphot)
		cphot=0.0_dp
		do i=1,nphotpix
		  cphot=cphot+w(photpixels(i))
		enddo
	
		!compute median(S/N)
		medsnr=0.0_dp
		if (nov > 0) then 
		  call rmedian(obs,e_obs,nlambda1,medsnr)
		  if (status == 0) write(*,*)'median snr =',medsnr
		endif		

	        if (nov > 0) then
	            if (covprint .eq. 1) then 
			write(3,'(1x,a150,100(1x,ES12.4))')fname,opf,ospf,cphot,   &
        		medsnr,lchi,ocov
	            else 
			write(3,'(1x,a150,100(1x,F9.3))')fname,opf,ospf,cphot,   &
        		medsnr,lchi
	            endif
                endif 
                !$omp end critical

	endif !if on status > -10 and k>= nruns


	enddo	!loop on nruns

	if (j >= only_object(2)) then
		status=10
		call quit(status)
	endif	

	!reset status for next object, or exit
	status=0

	if (trkout /= 0) then
		close(11)
		if (trkout < 0) close(12)
	endif

	call ellapsed_time(etime)
	
	
enddo 
!$omp end do
!$omp end parallel

!close units
close(1)
close(2)
close(3)
close(4)
if (snr == -1._dp) close(5)
if((nfilter > 1 .or. cont > 0) .and. sffile.gt.' ') close(7)
if (lsf > 10) close(9)
if (f_access == 1) close(10)


!sort output files when nthreads>1
if (nthreads > 1) call fsort()

if ( allocated( pf) ) deallocate(pf)
if ( allocated( pf0) ) deallocate(pf0)
if ( allocated( spf) ) deallocate(spf)
if ( allocated( pt) ) deallocate(pt)
if ( allocated( obs) ) deallocate(obs)
if ( allocated( e_obs) ) deallocate(e_obs)
if ( allocated( w) ) deallocate(w)
if ( allocated( fit) ) deallocate(fit)
if ( allocated( sfit) ) deallocate(sfit)
if ( allocated( cov) ) deallocate(cov)
if ( allocated( ocov) ) deallocate(ocov)
if ( allocated( bestpf) ) deallocate(bestpf)
if ( allocated( bestspf) ) deallocate(bestspf)
if ( allocated( bestcov) ) deallocate(bestcov)
if ( allocated( opf) ) deallocate(opf)
if ( allocated( ospf) ) deallocate(ospf)
if ( allocated( lambda_obs) ) deallocate(lambda_obs)
if ( allocated( obs_in) ) deallocate(obs_in)
if ( allocated( e_obs_in) ) deallocate(e_obs_in)
if ( allocated( obspca) ) deallocate(obspca)
if ( allocated( e_obspca) ) deallocate(e_obspca)
if ( allocated( lsfcof1) ) deallocate(lsfcof1)
if ( allocated( lsfarr1) ) deallocate(lsfarr1)
if ( allocated( lambda_syn) ) deallocate(lambda_syn)
if ( allocated( lsfcof) ) deallocate(lsfcof)
if ( allocated( lsfarr) ) deallocate(lsfarr)
if (  allocated( probe) ) deallocate(probe)
if ( allocated( ee) ) deallocate(ee)
if ( allocated( aa) ) deallocate(aa)
if ( allocated( uu) ) deallocate(uu)
if ( allocated( imap) ) deallocate(imap)
if (  allocated( wref) ) deallocate(wref)
if ( allocated( waveline) ) deallocate(waveline)


end do

call quit(status)

end program ferre

