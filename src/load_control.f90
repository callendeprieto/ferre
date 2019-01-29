
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine	load_control

!reading control namelist
!this routine is to be run serial -- no openmp critical protection provided
		
use share
implicit none

!locals
integer :: i
		
namelist / lista / ndim,nov,indv
namelist / lista / ntie,indtie,typetie,ttie0,ttie
namelist / lista / indini,nobj,nlambda,synthfile,fixfile,filterfile
namelist / lista / pfile,ffile,erfile,wfile
namelist / lista / opfile,offile,sffile,lsffile
namelist / lista / f_format,f_access,fformat
namelist / lista / snr,only_object,ycutoff,wphot,balance
namelist / lista / optimize,impact,mforce,chiout,trkout,nfilter,init
namelist / lista / nruns,errbar,covprint,indi,winter,twinter
namelist / lista / inter,mono,algor,scope,stopcr,simp,nthreads
namelist / lista / pcaproject,pcachi,lsf,nlsf
namelist / lista / cont,ncont,obscont,rejectcont
namelist / lista / ext_chain_filename,ext_gr_filename
namelist / lista / chain_num, gen_num, burnin_limit


indini(1:maxndim)=-10
indi(1:maxndim)=-1
open(1,file='input.nml',delim='apostrophe',recl=siobuffer)
read(1,nml=lista)
close(1)

!basic checks of input
!maxndim>=ndim>=1
if (ndim > maxndim .or. ndim < 1) then
	write(*,*) 'ERROR in load_control'
	write(*,*) 'ndim = ',ndim,' is > maxndim =',maxndim,' or'
	write(*,*) 'ndim = ',ndim,' is < 1'
	stop
endif

!initialize indi with default values when not set in the control file
!indi=[ndim,ndim-1,ndim-2 ...,1,...]
if (indi(1) < 0) then 
  do i=1,ndim
    indi(i)=ndim-i+1
  enddo
endif

!further checks
!ndim>nov>=1
if (nov > ndim .or. nov < 1) then
	if (nov == 0) then 
	  write(*,*) 'load_control: WARNING'
	  write(*,*) 'nov=0, using Pure Interpolation mode'
	  nruns=1
	else
	  write(*,*) 'load_control: ERROR'
	  write(*,*) 'nov = ',nov,' is > ndim =',ndim,' or'
	  write(*,*) 'nov = ',nov,' is < 1 and not exactly 0'
	stop
	endif
endif
!ndim>=max(indv(1:nov)), min(indv(1:nov))>=1
if (maxval(indv(1:nov)) > ndim .or. minval(indv(1:nov)) < 1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'max(indv) = ',maxval(indv),' is > ndim =',ndim,' or'
	write(*,*) 'min(indv) = ',minval(indv),' is < 1'
	stop
endif
!1>f_format>0
if (f_format > 1 .or. f_format < 0) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'f_format = ',f_format,' can only be 0 or 1'
	stop
endif
!1>f_access>0
if (f_access > 1 .or. f_access < 0) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'f_access = ',f_access,' can only be 0 or 1'
	stop
endif
!1>fformat>0
if (fformat > 1 .or. fformat < 0) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'fformat = ',fformat,' can only be 0 or 1'
	stop
endif
!cont <=3
if (cont  > 3) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'cont = ',cont,' must be <=0, 1, 2 or 3'
	stop
endif
if (cont /= 0 .and. cont /= 1 .and. cont /= 2 .and. cont /=3) then
        write(*,*) 'load_control: ERROR'
        write(*,*) 'cont = ',cont,' must be 0,1,2 or 3'
        stop
endif
!ncont >=0
if (ncont < 0) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'ncont = ',ncont,' must be >=0'
	stop
endif
if (rejectcont <= 0.) then 
	write(*,*) 'load_control: ERROR'
	write(*,*) 'rejectcont = ',rejectcont,' must be >0'
	stop
endif
!mforce >=0 and <=2
if (mforce < 0 .or. mforce > 2) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'mforce = ',mforce,' must be 0, 1 or 2'
	stop
endif
!nfilter>=0
if (nfilter < 0) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'filter = ',nfilter,' is < 0'
	stop
endif

!default is start searching at grid center
!for nruns=1 (default) or at random for nruns>1
!but user may customize this behavior using inini(1:nov)
if (maxval(indini(1:nov)) == -10) then
	indini(1:nov)=1
	if (abs(nruns)>1) indini(1:nov)=0
endif
if (minval(indini(1:nov)) < -1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'indini is not defined for all variable  parameters'
	write(*,*) '(it contains entries < -1)'
	stop
endif
if (product(abs(indini(1:nov))) < nruns .and. product(abs(indini(1:nov))) > 0) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'nruns [=',nruns,']  > product(abs(indini(1:nov))) [=',product(abs(indini(1:nov))),']' 
	write(*,*) 'there is not enough starting points! '
	stop
endif 
if (product(abs(indini(1:nov))) .ne. nruns .and. product(abs(indini(1:nov))) > 0) then
	write(*,*) 'load_control: WARNING'
	write(*,*) 'nruns [=',nruns,']  < product(abs(indini(1:nov))) [=',product(abs(indini(1:nov))),']'
	write(*,*) 'not all the defined starting points will be used!'
endif

!3>=errbar>=0
if (errbar > 3 .or. (errbar < 0 .and. errbar /= -2)) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'errbar = ',errbar,' should be 0, 1, 2, or 3'
	stop
endif
if (errbar == 2 .and. abs(nruns)<2) then 
	write(*,*) 'load_control:  ERROR'
	write(*,*) 'errbar = ',errbar,'(nrunsigma) requires abs(nruns)>1'
	stop
endif
if (errbar == 2 .and. nruns < 2) then
        write(*,*) 'load_control:  WARN'
        write(*,*) 'To get error bars using errbar=-2 (nrunsigma)  you should set nruns>1'
endif
!winter can only be 0, 1 or 2
if (winter < 0 .or. winter > 2) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'winter = ',winter,' must be 0, 1 or 2'
	stop
endif
!winter > 0 needs an wfile
if (winter > 0 .and. wfile .le. '') then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'winter = ',winter,' is > 0'
	write(*,*) 'but no wfile is defined'
	stop
endif
!twinter can only be 0 or 1
if (twinter < 0 .or. twinter > 1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'twinter = ',twinter,' must be 0 or 1'
	stop
endif
!inter>=-1
if (inter < -1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'inter = ',inter,' is < -1'
	stop
endif
!mono == 0 or 1
if (mono /= 0 .and. mono /= 1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'mono = ',mono,' must be 0 or 1'
	stop
endif
if (mono == 1 .and. inter /= 3) then
        write(*,*) 'load_control: WARNING'
        write(*,*) 'mono = ',mono,' is only active for inter=3!'
endif
!-2<=algor<=4
if (algor < -1 .or. algor > 5) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'algor = ',algor,' has an ilegal value (should be -1,0,1,2,3 or 4)' 
	stop
endif
!nov=0 or nov>1 for algor=4
if (algor == 4 .and. nov==1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'This algorithm (Truncated Newton) cannot work with only 1 dimension (nov=1)'
	stop
endif
!pcachi and pcaproject
if (pcachi == 1 .and. pcaproject == 0) then
	write(*,*) 'load_control: WARNING'
	write(*,*) 'pcachi = 1 but pcaproject = 0, so pre-compressed observations are expected'
endif
if(pcachi == 0. .and. pcaproject == 1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'pcaproject = ',pcaproject,' must be set to 0 when pcachi = 0'
	stop
endif
!lsf
if (minval(abs(lsf-(/0,1,2,3,4,11,12,13,14/))) > 0) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'lsf = ',lsf,' is not one of the following: 0,1,2,3,4,11,12,13,14'
	stop
endif
if (lsf > 0 .and. lsffile .eq. '') then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'lsffile should be given when lsf > 0'
	stop
endif
if (lsf == 0 .and. lsffile .ne. '') then
	write(*,*) 'load_control: WARNING'
	write(*,*) 'lsffile is given but lsf = 0, so no convolution will be performed'
endif
if (minval(abs(lsf-(/1,2,11,12/))) == 0 .and. nlsf < 3) then 
	write(*,*) 'load_control: ERROR'
	write(*,*) 'lsf=',lsf,'  but nlsf =',nlsf,' is < 3'
	stop
endif
!mandatory input files
if (pfile .eq. '' .or. offile .eq. '') then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'pfile and offile should be given inputs'
	stop
endif
if (nov .ne. 0 .and. (offile .eq. '' .or. ffile .eq. '')) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'offile and ffile should be given inputs'
	stop
endif
if (nov .ne. 0 .and. snr < 0.0_dp .and. erfile .eq. '') then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'erfile or snr should be provided as input'
	stop
endif
if ((only_object(1) > 0 .or. only_object(2) < 10**lmaxnobj) .and. nthreads /= 1) then
	write(*,*) 'load_control: ERROR'
	write(*,*) 'only_object cannot be used in combination with openmp'
        write(*,*)'nthreads=',nthreads
	stop
endif
if (burnin_limit < 0) burnin_limit = gen_num/2
		
write(*,*)'load_control -> ndim=',ndim
write(*,*)'load_control -> nov =',nov
write(*,*)'load_control -> snr =',snr
write(*,*)'load_control -> nruns=',nruns
write(*,*)'load_control -> lsf=',lsf
write(*,*)'load_control -> nthreads=',nthreads

if (algor == 5) then
 write(*,*)'load_control -> gen_num=',gen_num
 write(*,*)'load_control -> burnin_limit=',burnin_limit
endif

		
end subroutine load_control
