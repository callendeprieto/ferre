module	share


implicit none
save

!all variables in this module are to be shared by omp threads

!parameters
!integer, parameter 	:: dp = selected_real_kind(14, 60)	!precision
integer, parameter	:: dp = selected_real_kind(6)
integer, parameter	:: maxndim=20			!limit to #dim
integer, parameter	:: maxsynth=1000 		!limit to #synth
integer, parameter  	:: maxnpca=1000			!limit to #npca sections
integer, parameter 	:: lmaxnobj=9			!limit to #objs
integer, parameter 	:: long = selected_int_kind(9)	!long integers
integer, parameter 	:: longenough = selected_int_kind(lmaxnobj)	!is 10**lmaxnobj
integer, parameter      :: flen=300 ! chars in strings for paths/files
real(dp), parameter :: lambdatol = 1.e-3_dp	!accepted wavelength error 
real(dp), parameter :: pi=3.1415926535897932384626433832795_dp

character(len=12)    	:: ver = 'v4.8.4'  !version


!params to read or built from synthfile 
integer			:: npca(maxnpca)	!# of input pixels to npca sections 
integer			:: nelnpca = 0, totalnpca = 0 !# of pca sections,sum(npca)
integer			:: nvar = 0 !# pca components per section 
real(dp),allocatable 	:: meanspca(:),vpca(:),wpca(:,:), ff(:,:) !pca arrays 
real(dp)        :: constant=0.0_dp	    ! a constant added to all data	
integer			:: n_p(maxndim)				!grid's size
integer		    :: ntot	!ntot=Prod_1^ndim(n_p)=n_p(1)*n_p(2)...
integer,allocatable     :: ntimes(:) !ntimes(i)=Prod_(i+1)^ndim(n_p(i))
integer			:: npix	= 0				!# of frequencies in library
real(dp)		:: llimits(maxndim),steps(maxndim)	!physical pars
integer			:: nphotpix=0				!# of phot. pix
integer,allocatable	    :: photpixels(:)	!photometry
real(dp),allocatable    :: lambda_syn(:) 	!wavelength array for library
real(dp)		:: scalef=1.0_dp	        !(<f>)
integer 		:: scaled = 0				!scalef applied to read_f
integer			:: transposed = 0           !f array is (npix,ntot) for transposed=0
									 !        of (ntot,npix) for transpose=1
character(len=flen) :: file_data19 ='' ! name of the atomic .19 linelist
character(len=flen) :: file_data20 ='' ! name of the moleculer .20 linelist
integer        	   	:: siobuffer=2048       !bytes
integer         	:: liobuffer=20000      !reset in read_f
integer			:: xliobuffer=20000	!reset in ascii2bin
integer			:: rango	  !dec. exponent range for dp = range(f) -- set in read_f
real(dp)		:: minusculo  !smallest number >0 for dp = tiny(f)   -- set in read_f


!params to keep track of whats in the different synth modules used
!these are used to separate spectral ranges and photometry at smooth1/2
type headsynth
real(dp)        :: res			!max. resolution
integer(long)   :: pixbegin,pixend	!boundaries
real(dp)        :: lambda0              !wave(1)
real(dp)	:: lambda1		!wave(2)
integer         :: lws			!equidist in log10lambda
real(dp)	:: lambdamin		!min(wavelength) AA
real(dp)	:: lambdamax		!max(wavelegnth) AA
end type headsynth

integer         :: nsynth		!# of modules
type(headsynth)	:: hs(maxsynth)


! params to be read from control file 
integer			:: ndim	!dim of grid
integer		  	:: nov	!# of variable parameters (0 for pure interpolation)
integer		  	:: indv(maxndim)  !indices of the var pars.
integer			:: ntie = 1 !# of tied parameters 
integer			:: indtie(maxndim)!indices of the tied pars.
integer			:: typetie = 0    !0 linear ties with pars.
					  !1 linear ties with param. deltas
real(dp)                :: ttie0(maxndim)
real(dp)                :: ttie(maxndim,maxndim) !arrays with coeff. for ties
			 ! p(indtie(j))= ttie0(j)+sum(ttie(j,1:ndim)*p(1:ndim))
			 ! j=1,ntie
integer			:: indini(maxndim)!init type for var pars.
                         !<0 start at value in pfile 
			 !0 start at random
			 !1 start at grid center
			 !>1 start at the center of indini(j) equidistant cells 
		         !for variable j=1,...,nov
integer         	:: nlambda  = 0   !# of frequencies in the input spectra
integer         	:: nlambda1 = 0   !actual # of frequencies used in chi2 eval
integer(longenough)	:: nobj	= 10**lmaxnobj	!number of objects
						!nobj<= 0 program counts them
character(len=flen) 	:: synthfile(maxsynth)  !grid file(s)
character(len=flen) 	:: fixfile(maxsynth)	!file(s) for flux ratio corrections
character(len=flen)     :: filterfile=''!file with input reference weights
character(len=flen) 	:: pfile=''		!file with input pars
character(len=flen) 	:: ffile=''		!file with input fluxes
character(len=flen) 	:: erfile=''	!file with input flux errors
character(len=flen)     :: wfile=''     !file with input wavelengths
character(len=flen) 	:: opfile=''	!file with output pars
character(len=flen) 	:: offile=''	!file with output fluxes
character(len=flen)	:: sffile=''	!smoothed/normalized output fluxes
character(len=flen)	:: lsffile=''	!file with input lsf
integer         	:: f_format = 0     !format for synth files
						!0=ascii, 1=unformatted
integer         	:: f_access = 0     !access style for synth files
						!0=ram, 1=direct-access file
integer         	:: fformat = 0		!format for flux files 
						!0=flat, 1=iqpackage
real(dp)	  	:: snr	= -1.0_dp	!S/N (-1 => there is an erfile)
integer(longenough)	:: only_object(2) = (/0,10**lmaxnobj/)!objects to run
real(dp)		:: ycutoff = -huge(1.0_dp)	!exclude strong line core
real(dp)		:: wphot = -1.0_dp	!photom.weigth<0 ->npix/nphotpix
integer         	:: balance = 0      !use weights w(i)=1/obs(i)**2 to balance
integer			:: optimize = 0		!optimize weights
integer			:: impact = 0		!use impact factors to optimize
integer			:: cont = 0		!0=no normalization, 1=polynomial fit, 2=pem, 3=running mean
integer			:: ncont=0		!order/pieces/filter width -1 for cont=1/2/3 respectively
integer			:: obscont=1	    !if ncont>0 and obscont/=0 normalize both data and models
integer			:: mforce = 0       !force equal mean/median between obs and 
                                            !flux arrays in fun.f90 (1=force equal mean, 2=force median)
integer			:: chiout = 0	!output chi**2 surfaces
integer                 :: trkout = 0   !output tracking params 
					!(1=norm,2=phys,<0 to print flux res.)
integer			:: nfilter = 0	!boxcar filtering nfilter+1 wide
					!no filtering for nfilter<=1
integer			:: init = 1 !init sets starting values
				!  init=0 starts at pfile values
				!       1 starts following
				!	  instructions in indini 
				!         (default is center for
				!	   nruns=1 or random for 							   nruns>1)
				!       2 uses likelyhood solution
				!       3 starts at best-fitting pix
integer			:: nruns = 1		!runs per object 
						!abs(nruns)>1 => init=2(load_control)
						!nruns>0 keeps record of best fitting
						!nruns<0 keeps record of all fittings
integer			:: errbar = 0	!method for error determination
						!0=getsigma,1=covsigma,2=mcsigma
						!-2=nrunsigma	
						!3=pdfsigma 
                                                !NOTE:
                                                !when algor=5 the MCMC error bars are adopted regardless of errbar
integer			:: covprint = 0         !output the cov. matrix to opfile
integer			:: indi(maxndim)	!order for the interpolations	
integer			:: winter = 0		!interpolate in wavelength
										!0=no, 1=interp. obs, 2=interp. mdl
integer			:: twinter = 0		!type of wavelength interpolation [0/1]
						!0  -> linear
						!1  -> cubic splines
integer			:: inter = 1		!type of param. interpolation [0/4]
						! <=0 -> nearest neighbor
						! 1   -> linear interpolation
						! 2   -> quadratic Bezier interp.
						! 3   -> cubic Bezier interpolation
						! 4   -> cubic splines
integer			:: mono = 0		!force monotonic interpolation
					        !avoiding creation of extrema
						!when using inter=3

integer			:: algor = 1  !1  ->N-M Miller's implementation
					      !2  ->Boender-Timmer-Rinnoy Kan global algorithm 
					      !3  ->Powell's UOBYQA algorithm
					      !4  ->Nash's Truncated Newton algorithm
                                              !5  -> MCMC
					      !0  -> find best-fitting 'pixel' 
					      !-1 -> sum over param.space

integer			:: pcaproject=0     !1 project input data (errors too) 
						    !	for PCA grids
						    !0 input data is already projected
						    !   or does not need to be (pcachi=0)

integer			:: pcachi=0	    !1 eval chi2 in PCA space for PCA grids
						    !0 eval chi2 in flux space for PCA grids

						
integer			:: lsf = 0      !lsf for convolving the synth library
						!      IF Gaussian lsf, only FWHM is read from lsffile
						!      otherwise the whole profile is read
						! 0 -> no lsf convolution
						!the following options have 1 lsf for all objects
						! 1 -> lsf is 1D (not changing with lambda), one for all
						! 2 -> lsf is 2D (changing with lambda), one for all
						! 3 -> lsf is 1D and Gaussian 
						! 4 -> lsf is 2D (changing with lambda) and Gaussian
						!the following options have lsf changing for each object
						!11 -> lsf is 1D and particular for each object 
						!12 -> lsf is 2D and particular for each object
						!13 -> lsf is 1D Gaussian and particular for each object
						!14 -> lsf is 2D Gaussian and particular for each object
integer			:: mlsf = 1		!# of wavelengths for the lsf (1 or npix)
integer			:: nlsf = 1		!pixels for the lsf (at a given wavelength)
integer			:: dlsf = 1     !number of coefficients for analytical lsf
real(dp),allocatable:: lsfcof(:,:)  !holder for lsf coefficients (mlsf,dlsf)
real(dp),allocatable:: lsfarr(:,:)  !holder for a generic lsf (mlsf,nlsf)
									!note that this should be flipped in lambda
									!so ready for convolution 
						

!parameters for all optimization methods
integer, parameter 	:: iprint=-10    	!degree of verbose (quiet=-10)
real(dp)		:: scope = 0.45_dp	!search scope (normalized units)					
real(dp)		:: stopcr=1.e-4_dp  	!convergence criterion(1.e-7_dp)
integer, parameter	:: maxf=1000		!max # of objfun evaluations

!specific Nelder-Mead parameters from control file	
real(dp)        	:: simp=1.e-4_dp    !requested precision (1.e-6_dp)												

!other parameters for minim (Nelder-Mead minimization)
!shared between main and the caller routines (getmin) 
integer, parameter	:: nloop=8			!# of internal loops
integer, parameter	:: iquad=1			!switch on surface fitting
								  
!params to share between main, objfun, and getsigma/mcsigma
real(dp), allocatable	:: probe(:) 		! scales to get error bars

!params for MCMC
!first the ones kind of fixed
integer ( kind = 4 ),parameter    :: cr_num=3
real ( kind = 8 ),parameter       :: gr_threshold = 1.2
integer ( kind = 4 ),parameter    :: jumpstep = 5
integer ( kind = 4 ),parameter    :: pair_num = 3
integer ( kind = 4 ),parameter    :: printstep = 10
character ( len = 255 ),parameter :: restart_read_filename = ''
character ( len = 255 ),parameter :: restart_write_filename = ''

!now the MCMC parameters to be changed in the control file
character (len=flen)              :: ext_chain_filename=''  !e.g.'.chain000.dat'
character (len=flen)              :: ext_gr_filename = ''   ! e.g. '.gr'
integer ( kind = 4 )              :: chain_num = 10
integer ( kind = 4 )              :: gen_num = 500
integer ( kind = 4)               :: burnin_limit  = -1 !reset in load_control



!params to share between read_f and lin/qua/cub
real(dp), allocatable	:: f(:,:)		!synth grid

!params to share between ferre (main) and lin/minlocus 
integer,allocatable	:: ee(:,:)		!see getee
integer,allocatable	:: aa(:,:)		!see getaa
integer,allocatable	:: uu(:,:)		!see getuu
integer,allocatable     :: imap(:)


!params to share between read_f and lin/qua/cub/cova/getsigma
real(dp)	    :: badflux=-1000._dp	!invalid fluxes
character(len=30)   :: fmtformat !format string for the fmt output


!other omp-shared variables
real(dp), allocatable	:: wref(:)		!reference weights
real(dp), allocatable	:: waveline(:)	!iq-formatted wavelengths line

!params for omp 
integer         	:: nthreads = 1 !OMP_NUM_THREADS is used unless
					!specified in control file
					!0= uses all available

end module share
