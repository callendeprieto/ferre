MODULE trn

use share, only: dp,ndim,nlambda1,mlsf,nlsf,   simp
use fun

IMPLICIT NONE


REAL (dp), ALLOCATABLE, SAVE  :: gv(:), r(:), zk(:), v(:), sk(:), yk(:),  &
                                 diagb(:), sr(:), yr(:), hyr(:), hg(:),   &
                                 hyk(:), zsol(:), emat(:)

CONTAINS


SUBROUTINE sfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, n, x, f, g)

IMPLICIT NONE

real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)             ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray    
    
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f, g(:)
    
REAL (dp)  :: xh, fh, hinv, y(n)
INTEGER    :: i


hinv = 1.0_dp / simp
CALL objfun (w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,x,f)
y = x(1:n)
DO  i = 1,n
	xh   = x(i)
	y(i) = x(i) + simp
	CALL objfun (w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,y,fh)
 	g(i) = (fh - f) * hinv
 	y(i) = xh
END DO
RETURN    
    
END SUBROUTINE sfun


!%% TRUNCATED-NEWTON METHOD:  SUBROUTINES
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-14  Time: 16:13:02
 
!   FOR OTHER MACHINES, MODIFY ROUTINE MCHPR1 (MACHINE EPSILON)
!   WRITTEN BY:  STEPHEN G. NASH
!                OPERATIONS RESEARCH AND APPLIED STATISTICS DEPT.
!                GEORGE MASON UNIVERSITY
!                FAIRFAX, VA 22030
!******************************************************************

SUBROUTINE lmqn (w,pf,pf0,obs,lambda_obs,e_obs, mobs,lsfarr, & 
				ifail, n, x, f, g, msglvl, maxit, maxfun, eta, &
                 stepmx, accrcy, xtol)
                 
                 
real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)                     ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray
                 

! N.B. Arguments W & LW have been removed.

INTEGER, INTENT(OUT)       :: ifail
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(OUT)     :: f
REAL (dp), INTENT(OUT)     :: g(n)
INTEGER, INTENT(IN)        :: msglvl
INTEGER, INTENT(IN)        :: maxit
INTEGER, INTENT(IN)        :: maxfun
REAL (dp), INTENT(IN)      :: eta
REAL (dp), INTENT(IN)      :: stepmx
REAL (dp), INTENT(IN)      :: accrcy
REAL (dp), INTENT(IN)      :: xtol

! THIS ROUTINE IS A TRUNCATED-NEWTON METHOD.
! THE TRUNCATED-NEWTON METHOD IS PRECONDITIONED BY A LIMITED-MEMORY
! QUASI-NEWTON METHOD (THIS PRECONDITIONING STRATEGY IS DEVELOPED
! IN THIS ROUTINE) WITH A FURTHER DIAGONAL SCALING (SEE ROUTINE NDIA3).
! FOR FURTHER DETAILS ON THE PARAMETERS, SEE ROUTINE TN.

INTEGER   :: i, icycle, ipivot(1), ireset,  &
             modet, nfeval, nftotl, niter, nlincg, nm1, nmodif, numf, nwhy
REAL (dp) :: abstol, alpha, difnew, difold, epsmch, epsred, fkeep, fm,  &
             fnew, ftest, gnorm, gsk, gtg, gtpnew, oldf, oldgtp, one, pe,  &
             peps, pnorm, reltol, rteps, rtleps, rtol, rtolsq, small, spe, &
             tnytol, toleps, xnorm, yksk, yrsr, zero
LOGICAL   :: lreset, upd1


! THE FOLLOWING IMSL AND STANDARD FUNCTIONS ARE USED

! EXTERNAL sfun
! COMMON /subscr/ lgv,lz1,lzk,lv,lsk,lyk,ldiagb,lsr,lyr,  &
!                 loldg,lhg,lhyk,lpk,lemat,lwtest

! INITIALIZE PARAMETERS AND CONSTANTS

IF (msglvl >= -2) WRITE(*,800)
CALL setpar(n)
upd1 = .true.
ireset = 0
nfeval = 0
nmodif = 0
nlincg = 0
zero = 0.d0
one = 1.d0
nm1 = n - 1

! WITHIN THIS ROUTINE THE ARRAY W(LOLDG) IS SHARED BY W(LHYR)

! lhyr = loldg

! CHECK PARAMETERS AND SET CONSTANTS

CALL chkucp(maxfun,nwhy,n,alpha,epsmch, eta,peps,rteps,rtol,  &
            rtolsq,stepmx,ftest,xtol,xnorm,x,small,accrcy)
IF (nwhy < 0) GO TO 120
CALL setucr(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,small,nftotl,niter,n,f,fnew,fm,gtg,oldf,g,x)
IF (msglvl >= 1) WRITE(*,810) niter,nftotl,nlincg,fnew,gtg

! CHECK FOR SMALL GRADIENT AT THE STARTING POINT.

ftest = one + ABS(fnew)
IF (gtg < 1.d-4*epsmch*ftest*ftest) GO TO 90

! SET INITIAL VALUES TO OTHER PARAMETERS

icycle = nm1
toleps = rtol + rteps
rtleps = rtolsq + epsmch
gnorm  = SQRT(gtg)
difnew = zero
epsred = 5.0D-2
fkeep  = fnew

! SET THE DIAGONAL OF THE APPROXIMATE HESSIAN TO UNITY.

diagb = one

! ..................START OF MAIN ITERATIVE LOOP..........

! COMPUTE THE NEW SEARCH DIRECTION

modet = msglvl - 3
CALL modlnp(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, & 
			modet,zsol,gv,r,v, diagb,emat,x,g,  &
            zk, n,niter,maxit,nfeval,nmodif, nlincg,upd1,yksk,  &
            gsk,yrsr,lreset,.false.,ipivot, accrcy,gtpnew,gnorm,xnorm)

20 hyr = g(1:n)
pnorm = dnrm2(n,zsol,1)
oldf = fnew
oldgtp = gtpnew

! PREPARE TO COMPUTE THE STEP LENGTH

pe = pnorm + epsmch

! COMPUTE THE ABSOLUTE AND RELATIVE TOLERANCES FOR THE LINEAR SEARCH

reltol = rteps*(xnorm + one)/pe
abstol = - epsmch*ftest/(oldgtp - epsmch)

! COMPUTE THE SMALLEST ALLOWABLE SPACING BETWEEN POINTS IN THE LINEAR SEARCH

tnytol = epsmch*(xnorm + one)/pe
spe = stepmx/pe

! SET THE INITIAL STEP LENGTH.

alpha = step1(fnew,fm,oldgtp,spe)

! PERFORM THE LINEAR SEARCH

CALL linder(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,n,small,epsmch,reltol,abstol,tnytol,  &
            eta,zero,spe,zsol,oldgtp,x,fnew,alpha,g,numf, nwhy)

niter = niter + 1
nftotl = nftotl + numf
gtg = SUM( g(1:n)**2 )
IF (msglvl >= 1) WRITE(*,810) niter,nftotl,nlincg,fnew,gtg
IF (nwhy < 0) GO TO 120
IF (nwhy == 0 .OR. nwhy == 2) GO TO 30

! THE LINEAR SEARCH HAS FAILED TO FIND A LOWER POINT

nwhy = 3
GO TO 100

30 IF (nwhy <= 1) GO TO 40
CALL sfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, n,x,fnew,g)
nftotl = nftotl + 1

! TERMINATE IF MORE THAN MAXFUN EVALUTATIONS HAVE BEEN MADE

40 nwhy = 2
IF (nftotl > maxfun) GO TO 110
nwhy = 0

! SET UP PARAMETERS USED IN CONVERGENCE AND RESETTING TESTS

difold = difnew
difnew = oldf - fnew

! IF THIS IS THE FIRST ITERATION OF A NEW CYCLE, COMPUTE THE
! PERCENTAGE REDUCTION FACTOR FOR THE RESETTING TEST.

IF (icycle /= 1) GO TO 50
IF (difnew > 2.0D0 *difold) epsred = epsred + epsred
IF (difnew < 5.0D-1*difold) epsred = 5.0D-1*epsred

50 gnorm = SQRT(gtg)
ftest = one + ABS(fnew)
xnorm = dnrm2(n,x,1)

! TEST FOR CONVERGENCE

IF ((alpha*pnorm < toleps*(one + xnorm) .AND. ABS(difnew) < rtleps*ftest  &
    .AND. gtg < peps*ftest*ftest) .OR. gtg < 1.d-4*accrcy*ftest*ftest) GO TO 90

! COMPUTE THE CHANGE IN THE ITERATES AND THE CORRESPONDING CHANGE
! IN THE GRADIENTS

DO  i = 1,n
  yk(i) = g(i) - hyr(i)
  sk(i) = alpha*zsol(i)
END DO

! SET UP PARAMETERS USED IN UPDATING THE DIRECTION OF SEARCH.

yksk = DOT_PRODUCT( yk, sk )
lreset = .false.
IF (icycle == nm1 .OR. difnew < epsred*(fkeep-fnew)) lreset = .true.
IF (lreset) GO TO 70
yrsr = DOT_PRODUCT( yr, sr )
IF (yrsr <= zero) lreset = .true.

70 upd1 = .false.

!      COMPUTE THE NEW SEARCH DIRECTION

modet = msglvl - 3
CALL modlnp(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, & 
            modet,zsol,gv,r,v, diagb,emat,x,g, &
            zk, n,niter,maxit,nfeval,nmodif, nlincg,upd1,yksk, &
            gsk,yrsr,lreset,.false.,ipivot, accrcy,gtpnew,gnorm,xnorm)
IF (lreset) GO TO 80

!      STORE THE ACCUMULATED CHANGE IN THE POINT AND GRADIENT AS AN
!      "AVERAGE" DIRECTION FOR PRECONDITIONING.

sr = sr + sk
yr = yr + yk
icycle = icycle + 1
GO TO 20

! RESET

80  ireset = ireset + 1

! INITIALIZE THE SUM OF ALL THE CHANGES IN X.

sr = sk
yr = yk
fkeep = fnew
icycle = 1
GO TO 20

! ...............END OF MAIN ITERATION.......................

90 ifail = 0
f = fnew
RETURN

100 oldf = fnew

! LOCAL SEARCH HERE COULD BE INSTALLED HERE

110 f = oldf

! SET IFAIL

120 ifail = nwhy
RETURN

800 FORMAT(//' NIT   NF   CG         F', t46, 'GTG'//)
810 FORMAT(' ', i3, ' ', i4, ' ', i4, ' ', g22.15, '  ', g15.8)
END SUBROUTINE lmqn


SUBROUTINE lmqnbc (w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr, & 
				   ifail, n, x, f, g, low, up, ipivot, msglvl, &
                   maxit, maxfun, eta, stepmx, accrcy, xtol)

real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)                     ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray    
	
! N.B. Arguments W & LW have been removed.

INTEGER, INTENT(OUT)       :: ifail
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN OUT)  :: f
REAL (dp), INTENT(IN OUT)  :: g(n)
REAL (dp), INTENT(IN OUT)  :: low(n)
REAL (dp), INTENT(IN OUT)  :: up(n)
INTEGER, INTENT(IN OUT)    :: ipivot(n)
INTEGER, INTENT(IN)        :: msglvl
INTEGER, INTENT(IN OUT)    :: maxit
INTEGER, INTENT(IN)        :: maxfun
REAL (dp), INTENT(IN OUT)  :: eta
REAL (dp), INTENT(IN OUT)  :: stepmx
REAL (dp), INTENT(IN OUT)  :: accrcy
REAL (dp), INTENT(IN OUT)  :: xtol


! THIS ROUTINE IS A BOUNDS-CONSTRAINED TRUNCATED-NEWTON METHOD.
! THE TRUNCATED-NEWTON METHOD IS PRECONDITIONED BY A LIMITED-MEMORY
! QUASI-NEWTON METHOD (THIS PRECONDITIONING STRATEGY IS DEVELOPED
! IN THIS ROUTINE) WITH A FURTHER DIAGONAL SCALING (SEE ROUTINE NDIA3).
! FOR FURTHER DETAILS ON THE PARAMETERS, SEE ROUTINE TNBC.

INTEGER   :: i, icycle, nftotl, niter, nm1, numf, nwhy
REAL (dp) :: abstol, alpha, difnew, difold, epsmch, epsred, fkeep, flast, &
             fm, fnew, ftest, gnorm, gsk, gtg, gtpnew, oldf, oldgtp, one, &
             pe, peps, pnorm, reltol, rteps, rtleps, rtol, rtolsq, small, &
             spe, tnytol, toleps, xnorm, yksk, yrsr, zero
LOGICAL   :: conv, lreset, upd1, newcon

INTEGER    :: ier, ireset, modet, nfeval, nlincg, nmodif
! EXTERNAL sfun
! COMMON/subscr/ lgv, lz1, lzk, lv, lsk, lyk, ldiagb, lsr, lyr,  &
!                loldg, lhg, lhyk, lpk, lemat, lwtest

! CHECK THAT INITIAL X IS FEASIBLE AND THAT THE BOUNDS ARE CONSISTENT

CALL crash(n,x,ipivot,low,up,ier)
IF (ier /= 0) WRITE(*,800)
IF (ier /= 0) RETURN
IF (msglvl >= 1) WRITE(*,810)

! INITIALIZE VARIABLES

CALL setpar(n)
upd1 = .true.
ireset = 0
nfeval = 0
nmodif = 0
conv = .false.
zero = 0.d0
one = 1.d0
nm1 = n - 1
nlincg = 0

! WITHIN THIS ROUTINE THE ARRAY W(LOLDG) IS SHARED BY W(LHYR)

! lhyr = loldg

! CHECK PARAMETERS AND SET CONSTANTS

CALL chkucp(maxfun,nwhy,n,alpha,epsmch,eta,peps,rteps,rtol,rtolsq, &
            stepmx,ftest,xtol,xnorm,x,small,accrcy)
IF (nwhy < 0) GO TO 160
CALL setucr(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,small,nftotl,niter,n,f,fnew,fm,gtg,oldf,g,x)
flast = fnew

! TEST THE LAGRANGE MULTIPLIERS TO SEE IF THEY ARE NON-NEGATIVE.
! BECAUSE THE CONSTRAINTS ARE ONLY LOWER BOUNDS, THE COMPONENTS
! OF THE GRADIENT CORRESPONDING TO THE ACTIVE CONSTRAINTS ARE THE
! LAGRANGE MULTIPLIERS.  AFTERWORDS, THE PROJECTED GRADIENT IS FORMED.

DO  i = 1,n
  IF (ipivot(i) == 2) CYCLE
  IF (-ipivot(i)*g(i) >= 0.d0) CYCLE
  ipivot(i) = 0
END DO
CALL ztime(n,g,ipivot)
gtg = SUM( g(1:n)**2 )
IF (msglvl >= 1) CALL monit(n, x, fnew, g, niter, nftotl, nfeval, ipivot)

! CHECK IF THE INITIAL POINT IS A LOCAL MINIMUM.

ftest = one + ABS(fnew)
IF (gtg < 1.d-4*epsmch*ftest*ftest) GO TO 130

! SET INITIAL VALUES TO OTHER PARAMETERS

icycle = nm1
toleps = rtol + rteps
rtleps = rtolsq + epsmch
gnorm  = SQRT(gtg)
difnew = zero
epsred = 5.0D-2
fkeep  = fnew

! SET THE DIAGONAL OF THE APPROXIMATE HESSIAN TO UNITY.

diagb = one

! ..................START OF MAIN ITERATIVE LOOP..........

! COMPUTE THE NEW SEARCH DIRECTION

modet = msglvl - 3
CALL modlnp(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, & 
            modet, zsol, gv, r, v, diagb, emat, x, g, &
            zk, n, niter, maxit, nfeval, nmodif, nlincg, upd1,  &
            yksk, gsk, yrsr, lreset, .true., ipivot, accrcy, gtpnew,  &
            gnorm, xnorm)

20 hyr = g(1:n)
pnorm = dnrm2(n,zsol,1)
oldf = fnew
oldgtp = gtpnew

! PREPARE TO COMPUTE THE STEP LENGTH

pe = pnorm + epsmch

! COMPUTE THE ABSOLUTE AND RELATIVE TOLERANCES FOR THE LINEAR SEARCH

reltol = rteps*(xnorm + one)/pe
abstol = - epsmch*ftest/(oldgtp - epsmch)

! COMPUTE THE SMALLEST ALLOWABLE SPACING BETWEEN POINTS IN THE LINEAR SEARCH

tnytol = epsmch*(xnorm + one)/pe
CALL stpmax(stepmx,pe,spe,n,x,zsol,ipivot,low,up)

! SET THE INITIAL STEP LENGTH.

alpha = step1(fnew,fm,oldgtp,spe)

! PERFORM THE LINEAR SEARCH

CALL linder(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,n,small,epsmch,reltol,abstol,tnytol,  &
            eta,zero,spe,zsol,oldgtp,x,fnew,alpha,g,numf,nwhy)
newcon = .false.
IF (ABS(alpha-spe) > 10*epsmch) GO TO 30
newcon = .true.
nwhy   = 0
CALL modz(n,x,zsol,ipivot,epsmch,low,up)
flast = fnew

30 IF (msglvl >= 3) WRITE(*,820) alpha,pnorm
niter = niter + 1
nftotl = nftotl + numf

! IF REQUIRED, PRINT THE DETAILS OF THIS ITERATION

IF (msglvl >= 1) CALL monit(n, x, fnew, g, niter, nftotl, nfeval, ipivot)
IF (nwhy < 0) GO TO 160
IF (nwhy == 0 .OR. nwhy == 2) GO TO 40

! THE LINEAR SEARCH HAS FAILED TO FIND A LOWER POINT

nwhy = 3
GO TO 140
40 IF (nwhy <= 1) GO TO 50
CALL sfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, n,x,fnew,g)
nftotl = nftotl + 1

! TERMINATE IF MORE THAN MAXFUN EVALUATIONS HAVE BEEN MADE

50 nwhy = 2
IF (nftotl > maxfun) GO TO 150
nwhy = 0

! SET UP PARAMETERS USED IN CONVERGENCE AND RESETTING TESTS

difold = difnew
difnew = oldf - fnew

! IF THIS IS THE FIRST ITERATION OF A NEW CYCLE, COMPUTE THE
! PERCENTAGE REDUCTION FACTOR FOR THE RESETTING TEST.

IF (icycle /= 1) GO TO 60
IF (difnew > 2.d0*difold) epsred = epsred + epsred
IF (difnew < 5.0D-1*difold) epsred = 5.0D-1*epsred

60 gv = g(1:n)
CALL ztime(n,gv,ipivot)
gtg = SUM( gv**2 )
gnorm = SQRT(gtg)
ftest = one + ABS(fnew)
xnorm = dnrm2(n,x,1)

! TEST FOR CONVERGENCE

CALL cnvtst(conv,alpha,pnorm,toleps,xnorm,difnew,rtleps, ftest,gtg,  &
            peps,epsmch,gtpnew,fnew,flast,g,ipivot,n,accrcy)
IF (conv) GO TO 130
CALL ztime(n, g, ipivot)

! COMPUTE THE CHANGE IN THE ITERATES AND THE CORRESPONDING CHANGE
! IN THE GRADIENTS

IF (newcon) GO TO 90
DO  i = 1,n
  yk(i) = g(i) - hyr(i)
  sk(i) = alpha*zsol(i)
END DO

! SET UP PARAMETERS USED IN UPDATING THE PRECONDITIONING STRATEGY.

yksk = DOT_PRODUCT( yk, sk )
lreset = .false.
IF (icycle == nm1 .OR. difnew < epsred*(fkeep-fnew)) lreset = .true.
IF (lreset) GO TO 80
yrsr = DOT_PRODUCT( yr, sr )
IF (yrsr <= zero) lreset = .true.

80 upd1 = .false.

!      COMPUTE THE NEW SEARCH DIRECTION

90 IF (upd1 .AND. msglvl >= 3) WRITE(*,830)
IF (newcon .AND. msglvl >= 3) WRITE(*,840)
modet = msglvl - 3
CALL modlnp(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, & 
            modet,zsol,gv,r,v, diagb,emat,x,g,  &
            zk, n,niter,maxit,nfeval,nmodif, nlincg,upd1,yksk,  &
            gsk,yrsr,lreset,.true.,ipivot, accrcy,gtpnew,gnorm,xnorm)
IF (newcon) GO TO 20
IF (lreset) GO TO 110

! COMPUTE THE ACCUMULATED STEP AND ITS CORRESPONDING GRADIENT DIFFERENCE.

sr = sr + sk
yr = yr + yk
icycle = icycle + 1
GO TO 20

! RESET

110 ireset = ireset + 1

! INITIALIZE THE SUM OF ALL THE CHANGES IN X.

sr = sk
yr = yk
fkeep = fnew
icycle = 1
GO TO 20

! ...............END OF MAIN ITERATION.......................

130 ifail = 0
f = fnew
RETURN

140 oldf = fnew

! LOCAL SEARCH COULD BE INSTALLED HERE

150 f = oldf
IF (msglvl >= 1) CALL monit(n,x, f,g,niter,nftotl,nfeval,ipivot)

! SET IFAIL

160 ifail = nwhy
RETURN

800 FORMAT(' THERE IS NO FEASIBLE POINT; TERMINATING ALGORITHM')
810 FORMAT(//'  NIT   NF   CG         F', t47, 'GTG',//)
820 FORMAT('        LINESEARCH RESULTS:  ALPHA,PNORM', 2(g12.4))
830 FORMAT(' UPD1 IS TRUE - TRIVIAL PRECONDITIONING')
840 FORMAT(' NEWCON IS TRUE - CONSTRAINT ADDED IN LINESEARCH')
END SUBROUTINE lmqnbc


SUBROUTINE monit(n, x, f, g, niter, nftotl, nfeval, ipivot)

! N.B. Argument IRESET has been removed.

! PRINT RESULTS OF CURRENT ITERATION

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(n)
REAL (dp), INTENT(IN)  :: f
REAL (dp), INTENT(IN)  :: g(n)
INTEGER, INTENT(IN)    :: niter
INTEGER, INTENT(IN)    :: nftotl
INTEGER, INTENT(IN)    :: nfeval
INTEGER, INTENT(IN)    :: ipivot(n)

REAL (dp)  :: gtg
INTEGER    :: i

gtg = 0.d0
DO  i = 1,n
  IF (ipivot(i) /= 0) CYCLE
  gtg = gtg + g(i)*g(i)
END DO
WRITE(*,800) niter, nftotl, nfeval, f, gtg
RETURN

800 FORMAT(' ', i4, ' ', i4, ' ', i4, ' ', g22.15, '  ', g15.8)
END SUBROUTINE monit


SUBROUTINE ztime(n, x, ipivot)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
INTEGER, INTENT(IN)        :: ipivot(n)

! THIS ROUTINE MULTIPLIES THE VECTOR X BY THE CONSTRAINT MATRIX Z

INTEGER  :: i

DO  i = 1,n
  IF (ipivot(i) /= 0) x(i) = 0.d0
END DO
RETURN
END SUBROUTINE ztime


SUBROUTINE stpmax(stepmx, pe, spe, n, x, p, ipivot, low, up)

! COMPUTE THE MAXIMUM ALLOWABLE STEP LENGTH

REAL (dp), INTENT(IN)   :: stepmx
REAL (dp), INTENT(IN)   :: pe
REAL (dp), INTENT(OUT)  :: spe
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(n)
REAL (dp), INTENT(IN)   :: p(n)
INTEGER, INTENT(IN)     :: ipivot(n)
REAL (dp), INTENT(IN)   :: low(n)
REAL (dp), INTENT(IN)   :: up(n)

REAL (dp)  :: t
INTEGER    :: i

spe = stepmx / pe
! SPE IS THE STANDARD (UNCONSTRAINED) MAX STEP
DO  i = 1,n
  IF (ipivot(i) /= 0) CYCLE
  IF (p(i) == 0.d0) CYCLE
  IF (p(i) > 0.d0) GO TO 5
  t = low(i) - x(i)
  IF (t > spe*p(i)) spe = t / p(i)
  CYCLE

  5 t = up(i) - x(i)
  IF (t < spe*p(i)) spe = t / p(i)
END DO
RETURN
END SUBROUTINE stpmax



SUBROUTINE modz(n, x, p, ipivot, epsmch, low, up)

! N.B. Arguments FLAST & FNEW have been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN)      :: p(n)
INTEGER, INTENT(IN OUT)    :: ipivot(n)
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN)      :: low(n)
REAL (dp), INTENT(IN)      :: up(n)

REAL (dp)  :: tol
INTEGER    :: i

! UPDATE THE CONSTRAINT MATRIX IF A NEW CONSTRAINT IS ENCOUNTERED

DO  i = 1,n
  IF (ipivot(i) /= 0) CYCLE
  IF (p(i) == 0.d0) CYCLE
  IF (p(i) > 0.d0) GO TO 5
  tol = 10 * epsmch * (ABS(low(i)) + 1.d0)
  IF (x(i)-low(i) > tol) CYCLE
  ipivot(i) = -1
  x(i) = low(i)
  CYCLE

  5 tol = 10 * epsmch * (ABS(up(i)) + 1.d0)
  IF (up(i)-x(i) > tol) CYCLE
  ipivot(i) = 1
  x(i) = up(i)
END DO
RETURN
END SUBROUTINE modz


SUBROUTINE cnvtst(conv, alpha, pnorm, toleps, xnorm, difnew, rtleps, ftest, &
                  gtg, peps, epsmch, gtpnew, fnew, flast, g, ipivot, n, accrcy)

LOGICAL, INTENT(OUT)       :: conv
REAL (dp), INTENT(IN OUT)  :: alpha
REAL (dp), INTENT(IN)      :: pnorm
REAL (dp), INTENT(IN)      :: toleps
REAL (dp), INTENT(IN OUT)  :: xnorm
REAL (dp), INTENT(IN OUT)  :: difnew
REAL (dp), INTENT(IN OUT)  :: rtleps
REAL (dp), INTENT(IN OUT)  :: ftest
REAL (dp), INTENT(IN OUT)  :: gtg
REAL (dp), INTENT(IN OUT)  :: peps
REAL (dp), INTENT(IN OUT)  :: epsmch
REAL (dp), INTENT(IN OUT)  :: gtpnew
REAL (dp), INTENT(IN)      :: fnew
REAL (dp), INTENT(IN OUT)  :: flast
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: g(n)
INTEGER, INTENT(IN OUT)    :: ipivot(n)
REAL (dp), INTENT(IN OUT)  :: accrcy

LOGICAL    :: ltest
REAL (dp)  :: one, cmax, t
INTEGER    :: i, imax

! TEST FOR CONVERGENCE

imax = 0
cmax = 0.d0
ltest = flast - fnew <= -5.d-1*gtpnew
DO  i = 1,n
  IF (ipivot(i) == 0 .OR. ipivot(i) == 2) CYCLE
  t = -ipivot(i)*g(i)
  IF (t >= 0.d0) CYCLE
  conv = .false.
  IF (ltest) CYCLE
  IF (cmax <= t) CYCLE
  cmax = t
  imax = i
END DO
IF (imax == 0) GO TO 15
ipivot(imax) = 0
flast = fnew
RETURN

15 conv = .false.
one = 1.d0
IF ((alpha*pnorm >= toleps*(one + xnorm) .OR. ABS(difnew) >= rtleps*ftest  &
    .OR. gtg >= peps*ftest*ftest) .AND. gtg >= 1.d-4*accrcy*ftest*ftest) RETURN
conv = .true.

! FOR DETAILS, SEE GILL, MURRAY, AND WRIGHT (1981, P. 308) AND
! FLETCHER (1981, P. 116).  THE MULTIPLIER TESTS (HERE, TESTING
! THE SIGN OF THE COMPONENTS OF THE GRADIENT) MAY STILL NEED TO
! MODIFIED TO INCORPORATE TOLERANCES FOR ZERO.

RETURN
END SUBROUTINE cnvtst



SUBROUTINE crash(n, x, ipivot, low, up, ier)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
INTEGER, INTENT(OUT)       :: ipivot(n)
REAL (dp), INTENT(IN)      :: low(n)
REAL (dp), INTENT(IN)      :: up(n)
INTEGER, INTENT(OUT)       :: ier

! THIS INITIALIZES THE CONSTRAINT INFORMATION, AND ENSURES THAT THE
! INITIAL POINT SATISFIES  LOW <= X <= UP.
! THE CONSTRAINTS ARE CHECKED FOR CONSISTENCY.

INTEGER  :: i

ier = 0
DO  i = 1,n
  IF (x(i) < low(i)) x(i) = low(i)
  IF (x(i) > up(i)) x(i) = up(i)
  ipivot(i) = 0
  IF (x(i) == low(i)) ipivot(i) = -1
  IF (x(i) == up(i)) ipivot(i) = 1
  IF (up(i) == low(i)) ipivot(i) = 2
  IF (low(i) > up(i)) ier = -i
END DO
RETURN
END SUBROUTINE crash


! THE VECTORS SK AND YK, ALTHOUGH NOT IN THE CALL,
! ARE USED (VIA THEIR POSITION IN W) BY THE ROUTINE MSOLVE.

SUBROUTINE modlnp(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, & 
                  modet, zsol, gv, r, v, diagb, emat, x, g, zk, n,  &
                  niter, maxit, nfeval, nmodif, nlincg, upd1, yksk, gsk,   &
                  yrsr, lreset, bounds, ipivot, accrcy, gtp, gnorm,  &
                  xnorm)

real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)                     ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray    


! N.B. Arguments W & LW have been removed.

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: modet
REAL (dp), INTENT(OUT)     :: zsol(n)
REAL (dp), INTENT(IN OUT)  :: gv(n)
REAL (dp), INTENT(OUT)     :: r(n)
REAL (dp), INTENT(OUT)     :: v(n)
REAL (dp), INTENT(IN OUT)  :: diagb(n)
REAL (dp), INTENT(IN OUT)  :: emat(n)
REAL (dp), INTENT(IN)      :: x(n)
REAL (dp), INTENT(IN OUT)  :: g(n)
REAL (dp), INTENT(IN OUT)  :: zk(n)
INTEGER, INTENT(IN OUT)    :: niter
INTEGER, INTENT(IN)        :: maxit
INTEGER, INTENT(IN OUT)    :: nfeval
INTEGER, INTENT(IN OUT)    :: nmodif
INTEGER, INTENT(IN OUT)    :: nlincg
LOGICAL, INTENT(IN OUT)    :: upd1
REAL (dp), INTENT(IN OUT)  :: yksk
REAL (dp), INTENT(IN OUT)  :: gsk
REAL (dp), INTENT(IN OUT)  :: yrsr
LOGICAL, INTENT(IN OUT)    :: lreset
LOGICAL, INTENT(IN)        :: bounds
INTEGER, INTENT(IN OUT)    :: ipivot(1)
REAL (dp), INTENT(IN)      :: accrcy
REAL (dp), INTENT(OUT)     :: gtp
REAL (dp), INTENT(IN)      :: gnorm
REAL (dp), INTENT(IN OUT)  :: xnorm



REAL (dp)  :: alpha,beta,delta, pr, qold, qnew, qtest, rhsnrm, rz, &
              rzold, tol, vgv
INTEGER    :: i, k
LOGICAL    :: first
! EXTERNAL sfun

! THIS ROUTINE PERFORMS A PRECONDITIONED CONJUGATE-GRADIENT
! ITERATION IN ORDER TO SOLVE THE NEWTON EQUATIONS FOR A SEARCH
! DIRECTION FOR A TRUNCATED-NEWTON ALGORITHM.  WHEN THE VALUE OF THE
! QUADRATIC MODEL IS SUFFICIENTLY REDUCED, THE ITERATION IS TERMINATED.

! PARAMETERS

! MODET       - INTEGER WHICH CONTROLS AMOUNT OF OUTPUT
! ZSOL        - COMPUTED SEARCH DIRECTION
! G           - CURRENT GRADIENT
! GV,GZ1,V    - SCRATCH VECTORS
! R           - RESIDUAL
! DIAGB,EMAT  - DIAGONAL PRECONDITONING MATRIX
! NITER       - NONLINEAR ITERATION #
! FEVAL       - VALUE OF QUADRATIC FUNCTION

! *************************************************************
! INITIALIZATION
! *************************************************************

! GENERAL INITIALIZATION

IF (modet > 0) WRITE(*,800)
IF (maxit == 0) RETURN
first = .true.
rhsnrm = gnorm
tol = 1.d-12
qold = 0.d0

! INITIALIZATION FOR PRECONDITIONED CONJUGATE-GRADIENT ALGORITHM

CALL initpc(diagb,emat,n,modet, upd1,yksk,gsk,yrsr,lreset)
DO  i = 1,n
  r(i) = -g(i)
  v(i) = 0.d0
  zsol(i) = 0.d0
END DO

! ************************************************************
! MAIN ITERATION
! ************************************************************

DO  k = 1,maxit
  nlincg = nlincg + 1
  IF (modet > 1) WRITE(*,810) k
  
! CG ITERATION TO SOLVE SYSTEM OF EQUATIONS
  
  IF (bounds) CALL ztime(n,r,ipivot)
  CALL msolve(r,zk,n,upd1,yksk,gsk, yrsr,lreset,first)
  IF (bounds) CALL ztime(n,zk,ipivot)
  rz = DOT_PRODUCT( r(1:n), zk(1:n) )
  IF (rz/rhsnrm < tol) GO TO 80
  IF (k == 1) beta = 0.d0
  IF (k > 1) beta = rz/rzold
  v(1:n) = zk(1:n) + beta*v(1:n)
  IF (bounds) CALL ztime(n,v,ipivot)
  CALL gtims(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,v,gv,n,x,g,first,delta,accrcy,xnorm)
  IF (bounds) CALL ztime(n,gv,ipivot)
  nfeval = nfeval + 1
  vgv = DOT_PRODUCT( v(1:n), gv(1:n) )
  IF (vgv/rhsnrm < tol) GO TO 50
  CALL ndia3(n,emat,v,gv,r,vgv,modet)
  
! COMPUTE LINEAR STEP LENGTH
  
  alpha = rz / vgv
  IF (modet >= 1) WRITE(*,820) alpha
  
! COMPUTE CURRENT SOLUTION AND RELATED VECTORS
  
  zsol(1:n) = zsol(1:n) + alpha*v(1:n)
  r(1:n) = r(1:n) - alpha*gv(1:n)
  
! TEST FOR CONVERGENCE
  
  gtp = DOT_PRODUCT( zsol(1:n), g(1:n) )
  pr = DOT_PRODUCT( r(1:n), zsol(1:n) )
  qnew = 5.d-1 * (gtp + pr)
  qtest = k * (1.d0 - qold/qnew)
  IF (qtest < 0.d0) GO TO 70
  qold = qnew
  IF (qtest <= 5.d-1) GO TO 70
  
! PERFORM CAUTIONARY TEST
  
  IF (gtp > 0) GO TO 40
  rzold = rz
END DO

! TERMINATE ALGORITHM

k = k-1
GO TO 70

! TRUNCATE ALGORITHM IN CASE OF AN EMERGENCY

40 IF (modet >= -1) WRITE(*,830) k
zsol(1:n) = zsol(1:n) - alpha*v(1:n)
gtp = DOT_PRODUCT( zsol(1:n), g(1:n) )
GO TO 90

50 IF (modet > -2) WRITE(*,840)
IF (k > 1) GO TO 70

CALL msolve(g,zsol,n,upd1,yksk,gsk,yrsr,lreset,first)
CALL negvec(n,zsol)
IF (bounds) CALL ztime(n,zsol,ipivot)
gtp = DOT_PRODUCT( zsol(1:n), g(1:n) )

70 IF (modet >= -1) WRITE(*,850) k, gtp
GO TO 90

80 IF (modet >= -1) WRITE(*,860)
IF (k > 1) GO TO 70
zsol(1:n) = g(1:n)
CALL negvec(n,zsol)
IF (bounds) CALL ztime(n,zsol,ipivot)
gtp = DOT_PRODUCT( zsol(1:n), g(1:n) )
GO TO 70

! STORE (OR RESTORE) DIAGONAL PRECONDITIONING

90 diagb(1:n) = emat(1:n)
RETURN

800 FORMAT(' '//' ENTERING MODLNP')
810 FORMAT(' '//' ### ITERATION ',i2,' ###')
820 FORMAT(' ALPHA',g16.8)
830 FORMAT(' G(T)Z POSITIVE AT ITERATION ',i2, ' - TRUNCATING METHOD'/)
840 FORMAT(t12, 'HESSIAN NOT POSITIVE-DEFINITE')
850 FORMAT(/t9, 'MODLAN TRUNCATED AFTER ', i3, ' ITERATIONS',  &
           '  GTP = ', g14.6)
860 FORMAT(' PRECONDITIONING NOT POSITIVE-DEFINITE')
END SUBROUTINE modlnp


SUBROUTINE ndia3(n, e, v, gv, r, vgv, modet)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: e(n)
REAL (dp), INTENT(IN)      :: v(n)
REAL (dp), INTENT(IN)      :: gv(n)
REAL (dp), INTENT(IN)      :: r(n)
REAL (dp), INTENT(IN)      :: vgv
INTEGER, INTENT(IN)        :: modet

REAL (dp)  :: vr
INTEGER    :: i

! UPDATE THE PRECONDITIOING MATRIX BASED ON A DIAGONAL VERSION
! OF THE BFGS QUASI-NEWTON UPDATE.

vr = DOT_PRODUCT( v(1:n), r(1:n) )
DO  i = 1,n
  e(i) = e(i) - r(i)*r(i)/vr + gv(i)*gv(i)/vgv
  IF (e(i) > 1.d-6) CYCLE
  IF (modet > -2) WRITE(*,800) e(i)
  e(i) = 1.d0
END DO
RETURN

800 FORMAT(' *** EMAT NEGATIVE:  ',g16.8)
END SUBROUTINE ndia3


!      SERVICE ROUTINES FOR OPTIMIZATION

SUBROUTINE negvec(n,v)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: v(n)

! NEGATIVE OF THE VECTOR V

v(1:n) = -v(1:n)
RETURN
END SUBROUTINE negvec



SUBROUTINE lsout(iloc, itest, xmin, fmin, gmin, xw, fw, gw, u, a,  &
                 b, tol, eps, scxbd, xlamda)

INTEGER, INTENT(IN)    :: iloc
INTEGER, INTENT(IN)    :: itest
REAL (dp), INTENT(IN)  :: xmin
REAL (dp), INTENT(IN)  :: fmin
REAL (dp), INTENT(IN)  :: gmin
REAL (dp), INTENT(IN)  :: xw
REAL (dp), INTENT(IN)  :: fw
REAL (dp), INTENT(IN)  :: gw
REAL (dp), INTENT(IN)  :: u
REAL (dp), INTENT(IN)  :: a
REAL (dp), INTENT(IN)  :: b
REAL (dp), INTENT(IN)  :: tol
REAL (dp), INTENT(IN)  :: eps
REAL (dp), INTENT(IN)  :: scxbd
REAL (dp), INTENT(IN)  :: xlamda

! ERROR PRINTOUTS FOR GETPTC

REAL (dp) :: ya,yb,ybnd,yw,yu

yu = xmin + u
ya = a + xmin
yb = b + xmin
yw = xw + xmin
ybnd = scxbd + xmin
WRITE(*,800)
WRITE(*,810) tol,eps
WRITE(*,820) ya,yb
WRITE(*,830) ybnd
WRITE(*,840) yw,fw,gw
WRITE(*,850) xmin,fmin,gmin
WRITE(*,860) yu
WRITE(*,870) iloc,itest
RETURN

800 FORMAT(///' OUTPUT FROM LINEAR SEARCH')
810 FORMAT('  TOL AND EPS'/2g25.14)
820 FORMAT('  CURRENT UPPER AND LOWER BOUNDS'/2g25.14)
830 FORMAT('  STRICT UPPER BOUND'/g25.14)
840 FORMAT('  XW, FW, GW'/3g25.14)
850 FORMAT('  XMIN, FMIN, GMIN'/3g25.14)
860 FORMAT('  NEW ESTIMATE'/2g25.14)
870 FORMAT('  ILOC AND ITEST'/2I3)
END SUBROUTINE lsout


FUNCTION step1(fnew,fm,gtp,smax) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: fnew
REAL (dp), INTENT(IN)  :: fm
REAL (dp), INTENT(IN)  :: gtp
REAL (dp), INTENT(IN)  :: smax
REAL (dp)              :: fn_val

! ********************************************************
! STEP1 RETURNS THE LENGTH OF THE INITIAL STEP TO BE TAKEN ALONG THE
! VECTOR P IN THE NEXT LINEAR SEARCH.
! ********************************************************

REAL (dp) :: alpha, d, epsmch

epsmch = mchpr1()
d = ABS(fnew-fm)
alpha = 1.d0
IF (2.d0*d <= (-gtp) .AND. d >= epsmch) alpha = -2.d0*d/gtp
IF (alpha >= smax) alpha = smax
fn_val = alpha
RETURN
END FUNCTION step1



FUNCTION mchpr1() RESULT(fn_val)

! RETURNS THE VALUE OF EPSMCH, WHERE EPSMCH IS THE SMALLEST POSSIBLE
! REAL NUMBER SUCH THAT 1.0 + EPSMCH .GT. 1.0

REAL (dp)  :: fn_val

fn_val = EPSILON(0.0_dp)

RETURN
END FUNCTION mchpr1


SUBROUTINE chkucp(maxfun, nwhy, n, alpha, epsmch, eta, peps, rteps, &
                  rtol, rtolsq, stepmx, test, xtol, xnorm, x, small,   &
                  accrcy)

! N.B. Arguments LWTEST, LW & TINY has been removed.

INTEGER, INTENT(IN)     :: maxfun
INTEGER, INTENT(OUT)    :: nwhy
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: alpha
REAL (dp), INTENT(OUT)  :: epsmch
REAL (dp), INTENT(IN)   :: eta
REAL (dp), INTENT(OUT)  :: peps
REAL (dp), INTENT(OUT)  :: rteps
REAL (dp), INTENT(OUT)  :: rtol
REAL (dp), INTENT(OUT)  :: rtolsq
REAL (dp), INTENT(IN)   :: stepmx
REAL (dp), INTENT(OUT)  :: test
REAL (dp), INTENT(IN)   :: xtol
REAL (dp), INTENT(OUT)  :: xnorm
REAL (dp), INTENT(IN)   :: x(n)
REAL (dp), INTENT(OUT)  :: small
REAL (dp), INTENT(IN)   :: accrcy

! CHECKS PARAMETERS AND SETS CONSTANTS WHICH ARE COMMON TO BOTH
! DERIVATIVE AND NON-DERIVATIVE ALGORITHMS

epsmch = mchpr1()
small = epsmch*epsmch
nwhy = -1
rteps = SQRT(epsmch)
rtol = xtol
IF (ABS(rtol) < accrcy) rtol = 1.d1*rteps

! CHECK FOR ERRORS IN THE INPUT PARAMETERS

IF (n < 1 .OR. rtol < 0.d0 .OR. eta >= 1.d0 .OR.  &
    eta < 0.d0 .OR. stepmx < rtol .OR. maxfun < 1) RETURN
nwhy = 0

! SET CONSTANTS FOR LATER

rtolsq = rtol*rtol
peps = accrcy**0.6666_dp
xnorm = dnrm2(n,x,1)
alpha = 0.d0
test = 0.d0
RETURN
END SUBROUTINE chkucp


SUBROUTINE setucr(w,pf,pf0,obs,lambda_obs,e_obs, mobs,lsfarr, & 
				small, nftotl, niter, n, f, fnew, fm, gtg, oldf, g, x)
				
real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)                     ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray    				

REAL (dp), INTENT(IN)   :: small
INTEGER, INTENT(OUT)    :: nftotl
INTEGER, INTENT(OUT)    :: niter
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: f
REAL (dp), INTENT(OUT)  :: fnew
REAL (dp), INTENT(OUT)  :: fm
REAL (dp), INTENT(OUT)  :: gtg
REAL (dp), INTENT(OUT)  :: oldf
REAL (dp), INTENT(OUT)  :: g(n)
REAL (dp), INTENT(IN)   :: x(n)


! EXTERNAL         sfun

! CHECK INPUT PARAMETERS, COMPUTE THE INITIAL FUNCTION VALUE, SET
! CONSTANTS FOR THE SUBSEQUENT MINIMIZATION

fm = f

! COMPUTE THE INITIAL FUNCTION VALUE

CALL sfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, n,x,fnew,g)
nftotl = 1

! SET CONSTANTS FOR LATER

niter = 0
oldf = fnew
gtg = SUM( g(1:n)**2 )
RETURN
END SUBROUTINE setucr



SUBROUTINE gtims(w,pf,pf0,obs,lambda_obs,e_obs, mobs,lsfarr, & 
                 v, gv, n, x, g, first, delta, accrcy, xnorm)


real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)                     ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray    

! N.B. Arguments W & LW have been removed.

INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(IN)    :: v(n)
REAL (dp), INTENT(OUT)   :: gv(n)
REAL (dp), INTENT(IN)    :: x(n)
REAL (dp), INTENT(IN)    :: g(n)
LOGICAL, INTENT(IN OUT)  :: first
REAL (dp), INTENT(OUT)   :: delta
REAL (dp), INTENT(IN)    :: accrcy
REAL (dp), INTENT(IN)    :: xnorm


REAL (dp)  :: dinv
REAL (dp)  :: f
INTEGER    :: i

! EXTERNAL sfun
! COMMON/subscr/ lgv,lz1,lzk,lv,lsk,lyk,ldiagb,lsr,lyr,  &
!                lhyr,lhg,lhyk,lpk,lemat,lwtest

! THIS ROUTINE COMPUTES THE PRODUCT OF THE MATRIX G TIMES THE VECTOR
! V AND STORES THE RESULT IN THE VECTOR GV (FINITE-DIFFERENCE VERSION)

IF (.NOT. first) GO TO 20
delta = SQRT(accrcy)*(1.d0 + xnorm)
first = .false.

20 dinv = 1.d0/delta
hg = x(1:n) + delta*v(1:n)
CALL sfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, n,hg,f,gv)
DO  i = 1,n
  gv(i) = (gv(i) - g(i))*dinv
END DO
RETURN
END SUBROUTINE gtims



SUBROUTINE msolve(g, y, n, upd1, yksk, gsk, yrsr, lreset, first)

! N.B. Arguments W & LW have been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: g(n)
REAL (dp), INTENT(IN OUT)  :: y(n)
LOGICAL, INTENT(IN OUT)    :: upd1
REAL (dp), INTENT(IN OUT)  :: yksk
REAL (dp), INTENT(IN OUT)  :: gsk
REAL (dp), INTENT(IN OUT)  :: yrsr
LOGICAL, INTENT(IN OUT)    :: lreset
LOGICAL, INTENT(IN OUT)    :: first

! THIS ROUTINE SETS UPT THE ARRAYS FOR MSLV

! COMMON/subscr/ lgv,lz1,lzk,lv,lsk,lyk,ldiagb,lsr,lyr,  &
!                lhyr,lhg,lhyk,lpk,lemat,lwtest

CALL mslv(g,y,n,sk,yk,diagb,sr,yr,hyr,  &
          hg,hyk,upd1,yksk,gsk,yrsr,lreset,first)
RETURN
END SUBROUTINE msolve



SUBROUTINE mslv(g, y, n, sk, yk, diagb, sr, yr, hyr, hg, hyk,  &
                upd1, yksk, gsk, yrsr, lreset, first)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: g(n)
REAL (dp), INTENT(OUT)     :: y(n)
REAL (dp), INTENT(IN)      :: sk(n)
REAL (dp), INTENT(IN)      :: yk(n)
REAL (dp), INTENT(IN)      :: diagb(n)
REAL (dp), INTENT(IN)      :: sr(n)
REAL (dp), INTENT(IN)      :: yr(n)
REAL (dp), INTENT(OUT)     :: hyr(n)
REAL (dp), INTENT(OUT)     :: hg(n)
REAL (dp), INTENT(OUT)     :: hyk(n)
LOGICAL, INTENT(IN)        :: upd1
REAL (dp), INTENT(IN OUT)  :: yksk
REAL (dp), INTENT(OUT)     :: gsk
REAL (dp), INTENT(IN OUT)  :: yrsr
LOGICAL, INTENT(IN)        :: lreset
LOGICAL, INTENT(IN)        :: first

! THIS ROUTINE ACTS AS A PRECONDITIONING STEP FOR THE LINEAR CONJUGATE-
! GRADIENT ROUTINE.  IT IS ALSO THE METHOD OF COMPUTING THE SEARCH DIRECTION
! FROM THE GRADIENT FOR THE NON-LINEAR CONJUGATE-GRADIENT CODE.
! IT REPRESENTS A TWO-STEP SELF-SCALED BFGS FORMULA.

REAL (dp) :: rdiagb, ykhyk, ghyk, yksr, ykhyr, yrhyr, gsr, ghyr
REAL (dp) :: one = 1.0_dp
INTEGER   :: i

IF (upd1) GO TO 100
gsk = DOT_PRODUCT( g(1:n), sk(1:n) )
IF (lreset) GO TO 60

! COMPUTE HG AND HY WHERE H IS THE INVERSE OF THE DIAGONALS

DO  i = 1,n
  rdiagb = 1.0D0/diagb(i)
  hg(i) = g(i)*rdiagb
  IF (first) hyk(i) = yk(i)*rdiagb
  IF (first) hyr(i) = yr(i)*rdiagb
END DO
IF (first) THEN
  yksr  = DOT_PRODUCT( yk(1:n), sr(1:n) )
  ykhyr = DOT_PRODUCT( yk(1:n), hyr(1:n) )
  yrhyr = DOT_PRODUCT( yr(1:n), hyr(1:n) )
END IF
gsr  = DOT_PRODUCT( g(1:n), sr(1:n) )
ghyr = DOT_PRODUCT( g(1:n), hyr(1:n) )
CALL ssbfgs(n, one, sr, yr, hg, hyr, yrsr, yrhyr, gsr, ghyr, hg)
IF (first) CALL ssbfgs(n, one, sr, yr, hyk, hyr, yrsr, yrhyr, yksr, ykhyr, hyk)
ykhyk = DOT_PRODUCT( hyk(1:n), yk(1:n) )
ghyk  = DOT_PRODUCT( hyk(1:n), g(1:n) )
CALL ssbfgs(n, one, sk, yk, hg, hyk, yksk, ykhyk, gsk, ghyk, y)
RETURN

! COMPUTE GH AND HY WHERE H IS THE INVERSE OF THE DIAGONALS

60 DO  i = 1,n
  rdiagb = 1.d0/diagb(i)
  hg(i) = g(i)*rdiagb
  IF (first) hyk(i) = yk(i)*rdiagb
END DO
IF (first) ykhyk = DOT_PRODUCT( yk(1:n), hyk(1:n) )
ghyk = DOT_PRODUCT( g(1:n), hyk(1:n) )
CALL ssbfgs(n, one, sk, yk, hg, hyk, yksk, ykhyk, gsk, ghyk, y)
RETURN

100 y(1:n) = g(1:n) / diagb(1:n)
RETURN
END SUBROUTINE mslv


SUBROUTINE ssbfgs(n, gamma, sj, yj, hjv, hjyj, yjsj, yjhyj, vsj, vhyj, hjp1v)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: gamma
REAL (dp), INTENT(IN)   :: sj(n)
REAL (dp), INTENT(IN)   :: yj(n)
REAL (dp), INTENT(IN)   :: hjv(n)
REAL (dp), INTENT(IN)   :: hjyj(n)
REAL (dp), INTENT(IN)   :: yjsj
REAL (dp), INTENT(IN)   :: yjhyj
REAL (dp), INTENT(IN)   :: vsj
REAL (dp), INTENT(IN)   :: vhyj
REAL (dp), INTENT(OUT)  :: hjp1v(n)

! SELF-SCALED BFGS

INTEGER    :: i
REAL (dp)  :: beta,delta

delta = (1.d0 + gamma*yjhyj/yjsj)*vsj/yjsj - gamma*vhyj/yjsj
beta = -gamma*vsj/yjsj
DO  i = 1,n
  hjp1v(i) = gamma*hjv(i) + delta*sj(i) + beta*hjyj(i)
END DO
RETURN
END SUBROUTINE ssbfgs


! ROUTINES TO INITIALIZE PRECONDITIONER

SUBROUTINE initpc(diagb, emat, n, modet, upd1, yksk, gsk, yrsr, lreset)

! N.B. Arguments W & LW have been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: diagb(n)
REAL (dp), INTENT(IN OUT)  :: emat(n)
INTEGER, INTENT(IN)        :: modet
LOGICAL, INTENT(IN OUT)    :: upd1
REAL (dp), INTENT(IN OUT)  :: yksk
REAL (dp), INTENT(IN OUT)  :: gsk
REAL (dp), INTENT(IN OUT)  :: yrsr
LOGICAL, INTENT(IN OUT)    :: lreset

! COMMON/subscr/ lgv,lz1,lzk,lv,lsk,lyk,ldiagb,lsr,lyr,  &
!                lhyr,lhg,lhyk,lpk,lemat,lwtest

CALL initp3(diagb,emat,n,lreset,yksk,yrsr,hyk,  &
            sk,yk,sr,yr,modet,upd1)
RETURN
END SUBROUTINE initpc



SUBROUTINE initp3(diagb, emat, n, lreset, yksk, yrsr, bsk, sk, yk, sr, yr, &
                  modet, upd1)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: diagb(n)
REAL (dp), INTENT(OUT)  :: emat(n)
LOGICAL, INTENT(IN)     :: lreset
REAL (dp), INTENT(IN)   :: yksk
REAL (dp), INTENT(IN)   :: yrsr
REAL (dp), INTENT(OUT)  :: bsk(n)
REAL (dp), INTENT(IN)   :: sk(n)
REAL (dp), INTENT(IN)   :: yk(n)
REAL (dp), INTENT(IN)   :: sr(n)
REAL (dp), INTENT(IN)   :: yr(n)
INTEGER, INTENT(IN)     :: modet
LOGICAL, INTENT(IN)     :: upd1

REAL (dp)  :: cond, sds, srds, yrsk, td, d1, dn
INTEGER    :: i

IF (upd1) GO TO 90
IF (lreset) GO TO 60

bsk(1:n) = diagb(1:n)*sr(1:n)
sds  = DOT_PRODUCT( sr(1:n), bsk(1:n) )
srds = DOT_PRODUCT( sk(1:n), bsk(1:n) )
yrsk = DOT_PRODUCT( yr(1:n), sk(1:n) )
DO  i = 1,n
  td = diagb(i)
  bsk(i) = td*sk(i) - bsk(i)*srds/sds + yr(i)*yrsk/yrsr
  emat(i) = td - td*td*sr(i)*sr(i)/sds + yr(i)*yr(i)/yrsr
END DO
sds = DOT_PRODUCT( sk(1:n), bsk(1:n) )
DO  i = 1,n
  emat(i) = emat(i) - bsk(i)*bsk(i)/sds + yk(i)*yk(i)/yksk
END DO
GO TO 110

60 bsk(1:n) = diagb(1:n)*sk(1:n)
sds = DOT_PRODUCT( sk(1:n), bsk(1:n) )
DO  i = 1,n
  td = diagb(i)
  emat(i) = td - td*td*sk(i)*sk(i)/sds + yk(i)*yk(i)/yksk
END DO
GO TO 110

90 emat(1:n) = diagb(1:n)

110 IF (modet < 1) RETURN
d1 = emat(1)
dn = emat(1)
DO  i = 1,n
  IF (emat(i) < d1) d1 = emat(i)
  IF (emat(i) > dn) dn = emat(i)
END DO
cond = dn/d1
WRITE(*,800) d1,dn,cond
800 FORMAT(//t9, 'DMIN =', g12.4, '  DMAX =', g12.4, ' COND =', g12.4/)
RETURN
END SUBROUTINE initp3



SUBROUTINE setpar(n)

INTEGER, INTENT(IN)  :: n

! INTEGER  :: lsub(14)
! COMMON/subscr/ lsub,lwtest

! SET UP PARAMETERS FOR THE OPTIMIZATION ROUTINE

! DO  i = 1,14
!   lsub(i) = (i-1)*n + 1
! END DO
! lwtest = lsub(14) + n - 1

if (.not. allocated(gv)) then
ALLOCATE( gv(n), r(n), zk(n), v(n), sk(n), yk(n), diagb(n), sr(n), yr(n), &
          hyr(n), hg(n), hyk(n), zsol(n), emat(n) )
endif !CAP

RETURN
END SUBROUTINE setpar


!      LINE SEARCH ALGORITHMS OF GILL AND MURRAY

SUBROUTINE linder(w,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr, & 
				  n, small, epsmch, reltol, abstol, tnytol, eta,  &
                  sftbnd, xbnd, p, gtp, x, f, alpha, g, nftotl, iflag)

real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(inout) :: pf(ndim)                     ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray    
	
! N.B. Arguments W & LW have been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: small
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: reltol
REAL (dp), INTENT(IN OUT)  :: abstol
REAL (dp), INTENT(IN OUT)  :: tnytol
REAL (dp), INTENT(IN)      :: eta
REAL (dp), INTENT(IN OUT)  :: sftbnd
REAL (dp), INTENT(IN OUT)  :: xbnd
REAL (dp), INTENT(IN)      :: p(n)
REAL (dp), INTENT(IN)      :: gtp
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN OUT)  :: f
REAL (dp), INTENT(IN OUT)  :: alpha
REAL (dp), INTENT(IN OUT)  :: g(n)
INTEGER, INTENT(OUT)       :: nftotl
INTEGER, INTENT(OUT)       :: iflag

INTEGER    :: ientry, itest, lsprnt, numf, itcnt, nprnt
REAL (dp)  :: a, b,b1,big,e,factor,fmin,fpresn,fu, fw,gmin,gtest1,gtest2, &
              gu,gw,oldf,scxbnd,step, tol,u,xmin,xw,rmu,rtsmll,ualpha
LOGICAL    :: braktd
REAL (dp)  :: wx(n), wg(n)

!      THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE
!      CALLED WITHIN LINDER

! EXTERNAL sfun

!      ALLOCATE THE ADDRESSES FOR LOCAL WORKSPACE

lsprnt = 0
nprnt  = 10000
rtsmll = SQRT(small)
big = 1.d0/small
itcnt = 0

!      SET THE ESTIMATED RELATIVE PRECISION IN F(X).

fpresn = 10.d0*epsmch
numf = 0
u = alpha
fu = f
fmin = f
gu = gtp
rmu = 1.0D-4

!      FIRST ENTRY SETS UP THE INITIAL INTERVAL OF UNCERTAINTY.

ientry = 1

! TEST FOR TOO MANY ITERATIONS

10 itcnt = itcnt + 1
iflag = 1
IF (itcnt > 20) GO TO 50
iflag = 0
CALL getptc(big,rtsmll,reltol,abstol,tnytol, fpresn, eta, rmu,  &
            xbnd,u,fu,gu,xmin,fmin,gmin, xw,fw,gw,a,b,oldf,b1,scxbnd,  &
            e,step,factor, braktd,gtest1,gtest2,tol,ientry,itest)
!LSOUT
IF (lsprnt >= nprnt) CALL lsout(ientry,itest,xmin,fmin,gmin,  &
                                xw,fw,gw,u,a,b,tol,reltol,scxbnd,xbnd)

!      IF ITEST=1, THE ALGORITHM REQUIRES THE FUNCTION VALUE TO BE CALCULATED.

IF (itest /= 1) GO TO 30
ualpha = xmin + u
wx = x(1:n) + ualpha*p(1:n)
CALL sfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, n,wx,fu,wg)
numf = numf + 1
gu = DOT_PRODUCT( wg, p(1:n) )

!      THE GRADIENT VECTOR CORRESPONDING TO THE BEST POINT IS OVERWRITTEN
!      IF FU IS LESS THAN FMIN AND FU IS SUFFICIENTLY LOWER THAN F
!      AT THE ORIGIN.

IF (fu <= fmin .AND. fu <= oldf-ualpha*gtest1) g(1:n) = wg
GO TO 10

!      IF ITEST=2 OR 3 A LOWER POINT COULD NOT BE FOUND

30 nftotl = numf
iflag = 1
IF (itest /= 0) GO TO 50

!      IF ITEST=0 A SUCCESSFUL SEARCH HAS BEEN MADE

iflag = 0
f = fmin
alpha = xmin
x(1:n) = x(1:n) + alpha*p(1:n)

50 RETURN
END SUBROUTINE linder


SUBROUTINE getptc(big, rtsmll, reltol, abstol, tnytol, fpresn, eta,  &
                  rmu, xbnd, u, fu, gu, xmin, fmin, gmin, xw, fw, gw, a, b, &
                  oldf, b1, scxbnd, e, step, factor, braktd, gtest1, gtest2, &
                  tol, ientry, itest)

! N.B. Argument SMALL has been removed.

REAL (dp), INTENT(IN)      :: big
REAL (dp), INTENT(IN)      :: rtsmll
REAL (dp), INTENT(IN OUT)  :: reltol
REAL (dp), INTENT(IN OUT)  :: abstol
REAL (dp), INTENT(IN)      :: tnytol
REAL (dp), INTENT(IN)      :: fpresn
REAL (dp), INTENT(IN)      :: eta
REAL (dp), INTENT(IN)      :: rmu
REAL (dp), INTENT(IN)      :: xbnd
REAL (dp), INTENT(IN OUT)  :: u
REAL (dp), INTENT(IN OUT)  :: fu
REAL (dp), INTENT(IN OUT)  :: gu
REAL (dp), INTENT(OUT)     :: xmin
REAL (dp), INTENT(OUT)     :: fmin
REAL (dp), INTENT(OUT)     :: gmin
REAL (dp), INTENT(OUT)     :: xw
REAL (dp), INTENT(OUT)     :: fw
REAL (dp), INTENT(OUT)     :: gw
REAL (dp), INTENT(OUT)     :: a
REAL (dp), INTENT(OUT)     :: b
REAL (dp), INTENT(OUT)     :: oldf
REAL (dp), INTENT(OUT)     :: b1
REAL (dp), INTENT(OUT)     :: scxbnd
REAL (dp), INTENT(OUT)     :: e
REAL (dp), INTENT(OUT)     :: step
REAL (dp), INTENT(OUT)     :: factor
LOGICAL, INTENT(OUT)       :: braktd
REAL (dp), INTENT(OUT)     :: gtest1
REAL (dp), INTENT(OUT)     :: gtest2
REAL (dp), INTENT(OUT)     :: tol
INTEGER, INTENT(IN OUT)    :: ientry
INTEGER, INTENT(OUT)       :: itest

REAL (dp) :: denom

! ************************************************************
! GETPTC, AN ALGORITHM FOR FINDING A STEPLENGTH, CALLED REPEATEDLY BY
! ROUTINES WHICH REQUIRE A STEP LENGTH TO BE COMPUTED USING CUBIC
! INTERPOLATION.  THE PARAMETERS CONTAIN INFORMATION ABOUT THE INTERVAL
! IN WHICH A LOWER POINT IS TO BE FOUND AND FROM THIS GETPTC COMPUTES A
! POINT AT WHICH THE FUNCTION CAN BE EVALUATED BY THE CALLING PROGRAM.
! THE VALUE OF THE INTEGER PARAMETERS IENTRY DETERMINES THE PATH TAKEN
! THROUGH THE CODE.
! ************************************************************

LOGICAL    :: convrg
REAL (dp)  :: abgmin,abgw,absr,a1,chordm,chordu,  &
              d1,d2,p,q,r,s,scale,sumsq,twotol,xmidpt
REAL (dp)  :: zero, point1,half,one,three,five,eleven

! THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE CALLED
! WITHIN GETPTC

zero = 0.d0
point1 = 1.d-1
half = 5.d-1
one = 1.d0
three = 3.d0
five = 5.d0
eleven = 11.d0

!      BRANCH TO APPROPRIATE SECTION OF CODE DEPENDING ON THE VALUE OF IENTRY.

SELECT CASE ( ientry )
  CASE (    1)
    GO TO 10
  CASE (    2)
    GO TO 20
END SELECT

!      IENTRY=1
!      CHECK INPUT PARAMETERS

10 itest = 2
IF (u <= zero .OR. xbnd <= tnytol .OR. gu > zero) RETURN
itest = 1
IF (xbnd < abstol) abstol = xbnd
tol = abstol
twotol = tol + tol

! A AND B DEFINE THE INTERVAL OF UNCERTAINTY, X AND XW ARE POINTS
! WITH LOWEST AND SECOND LOWEST FUNCTION VALUES SO FAR OBTAINED.
! INITIALIZE A,SMIN,XW AT ORIGIN AND CORRESPONDING VALUES OF
! FUNCTION AND PROJECTION OF THE GRADIENT ALONG DIRECTION OF SEARCH
! AT VALUES FOR LATEST ESTIMATE AT MINIMUM.

a = zero
xw = zero
xmin = zero
oldf = fu
fmin = fu
fw = fu
gw = gu
gmin = gu
step = u
factor = five

!      THE MINIMUM HAS NOT YET BEEN BRACKETED.

braktd = .false.

! SET UP XBND AS A BOUND ON THE STEP TO BE TAKEN. (XBND IS NOT COMPUTED
! EXPLICITLY BUT SCXBND IS ITS SCALED VALUE.)  SET THE UPPER BOUND
! ON THE INTERVAL OF UNCERTAINTY INITIALLY TO XBND + TOL(XBND).

scxbnd = xbnd
b = scxbnd + reltol*ABS(scxbnd) + abstol
e = b + b
b1 = b

! COMPUTE THE CONSTANTS REQUIRED FOR THE TWO CONVERGENCE CRITERIA.

gtest1 = -rmu*gu
gtest2 = -eta*gu

! SET IENTRY TO INDICATE THAT THIS IS THE FIRST ITERATION

ientry = 2
GO TO 210

! IENTRY = 2

! UPDATE A,B,XW, AND XMIN

20 IF (fu > fmin) GO TO 60

! IF FUNCTION VALUE NOT INCREASED, NEW POINT BECOMES NEXT
! ORIGIN AND OTHER POINTS ARE SCALED ACCORDINGLY.

chordu = oldf - (xmin + u)*gtest1
IF (fu <= chordu) GO TO 30

! THE NEW FUNCTION VALUE DOES NOT SATISFY THE SUFFICIENT DECREASE
! CRITERION. PREPARE TO MOVE THE UPPER BOUND TO THIS POINT AND
! FORCE THE INTERPOLATION SCHEME TO EITHER BISECT THE INTERVAL OF
! UNCERTAINTY OR TAKE THE LINEAR INTERPOLATION STEP WHICH ESTIMATES
! THE ROOT OF F(ALPHA)=CHORD(ALPHA).

chordm = oldf - xmin*gtest1
gu = -gmin
denom = chordm-fmin
IF (ABS(denom) >= 1.d-15) GO TO 25
denom = 1.d-15
IF (chordm-fmin < 0.d0)  denom = -denom

25 IF (xmin /= zero) gu = gmin*(chordu-fu)/denom
fu = half*u*(gmin+gu) + fmin
IF (fu < fmin) fu = fmin
GO TO 60

30 fw = fmin
fmin = fu
gw = gmin
gmin = gu
xmin = xmin + u
a = a-u
b = b-u
xw = -u
scxbnd = scxbnd - u
IF (gu <= zero) GO TO 40
b = zero
braktd = .true.
GO TO 50

40 a = zero
50 tol = ABS(xmin)*reltol + abstol
GO TO 90

! IF FUNCTION VALUE INCREASED, ORIGIN REMAINS UNCHANGED
! BUT NEW POINT MAY NOW QUALIFY AS W.

60 IF (u < zero) GO TO 70
b = u
braktd = .true.
GO TO 80

70 a = u
80 xw = u
fw = fu
gw = gu
90 twotol = tol + tol
xmidpt = half*(a + b)

! CHECK TERMINATION CRITERIA

convrg = ABS(xmidpt) <= twotol - half*(b-a) .OR.  &
    ABS(gmin) <= gtest2 .AND. fmin < oldf .AND.  &
    (ABS(xmin - xbnd) > tol .OR. .NOT. braktd)
IF (.NOT. convrg) GO TO 100
itest = 0
IF (xmin /= zero) RETURN

! IF THE FUNCTION HAS NOT BEEN REDUCED, CHECK TO SEE THAT THE RELATIVE
! CHANGE IN F(X) IS CONSISTENT WITH THE ESTIMATE OF THE DELTA-
! UNIMODALITY CONSTANT, TOL.  IF THE CHANGE IN F(X) IS LARGER THAN
! EXPECTED, REDUCE THE VALUE OF TOL.

itest = 3
IF (ABS(oldf-fw) <= fpresn*(one + ABS(oldf))) RETURN
tol = point1*tol
IF (tol < tnytol) RETURN
reltol = point1*reltol
abstol = point1*abstol
twotol = point1*twotol

! CONTINUE WITH THE COMPUTATION OF A TRIAL STEP LENGTH

100 r = zero
q = zero
s = zero
IF (ABS(e) <= tol) GO TO 150

! FIT CUBIC THROUGH XMIN AND XW

r = three*(fmin-fw)/xw + gmin + gw
absr = ABS(r)
q = absr
IF (gw == zero .OR. gmin == zero) GO TO 140

! COMPUTE THE SQUARE ROOT OF (R*R - GMIN*GW) IN A WAY
! WHICH AVOIDS UNDERFLOW AND OVERFLOW.

abgw = ABS(gw)
abgmin = ABS(gmin)
s = SQRT(abgmin)*SQRT(abgw)
IF ((gw/abgw)*gmin > zero) GO TO 130

! COMPUTE THE SQUARE ROOT OF R*R + S*S.

sumsq = one
p = zero
IF (absr >= s) GO TO 110

! THERE IS A POSSIBILITY OF OVERFLOW.

IF (s > rtsmll) p = s*rtsmll
IF (absr >= p) sumsq = one +(absr/s)**2
scale = s
GO TO 120

! THERE IS A POSSIBILITY OF UNDERFLOW.

110 IF (absr > rtsmll) p = absr*rtsmll
IF (s >= p) sumsq = one + (s/absr)**2
scale = absr
120 sumsq = SQRT(sumsq)
q = big
IF (scale < big/sumsq) q = scale*sumsq
GO TO 140

! COMPUTE THE SQUARE ROOT OF R*R - S*S

130 q = SQRT(ABS(r+s))*SQRT(ABS(r-s))
IF (r >= s .OR. r <= (-s)) GO TO 140
r = zero
q = zero
GO TO 150

! COMPUTE THE MINIMUM OF FITTED CUBIC

140 IF (xw < zero) q = -q
s = xw*(gmin - r - q)
q = gw - gmin + q + q
IF (q > zero) s = -s
IF (q <= zero) q = -q
r = e
IF (b1 /= step .OR. braktd) e = step

! CONSTRUCT AN ARTIFICIAL BOUND ON THE ESTIMATED STEPLENGTH

150 a1 = a
b1 = b
step = xmidpt
IF (braktd) GO TO 160
step = -factor*xw
IF (step > scxbnd) step = scxbnd
IF (step /= scxbnd) factor = five*factor
GO TO 170

! IF THE MINIMUM IS BRACKETED BY 0 AND XW THE STEP MUST LIE WITHIN (A,B).

160 IF ((a /= zero .OR. xw >= zero) .AND. (b /= zero .OR.  &
    xw <= zero)) GO TO 180

! IF THE MINIMUM IS NOT BRACKETED BY 0 AND XW THE STEP MUST LIE WITHIN (A1,B1).

d1 = xw
d2 = a
IF (a == zero) d2 = b
! THIS LINE MIGHT BE
!     IF (A .EQ. ZERO) D2 = E
u = - d1/d2
step = five*d2*(point1 + one/u)/eleven
IF (u < one) step = half*d2*SQRT(u)
170 IF (step <= zero) a1 = step
IF (step > zero) b1 = step

! REJECT THE STEP OBTAINED BY INTERPOLATION IF IT LIES OUTSIDE THE
! REQUIRED INTERVAL OR IT IS GREATER THAN HALF THE STEP OBTAINED
! DURING THE LAST-BUT-ONE ITERATION.

180 IF (ABS(s) <= ABS(half*q*r) .OR. s <= q*a1 .OR. s >= q*b1) GO TO 200

! A CUBIC INTERPOLATION STEP

step = s/q

! THE FUNCTION MUST NOT BE EVALUTATED TOO CLOSE TO A OR B.

IF (step - a >= twotol .AND. b - step >= twotol) GO TO 210
IF (xmidpt > zero) GO TO 190
step = -tol
GO TO 210

190 step = tol
GO TO 210

200 e = b-a

! IF THE STEP IS TOO LARGE, REPLACE BY THE SCALED BOUND (SO AS TO
! COMPUTE THE NEW POINT ON THE BOUNDARY).

210 IF (step < scxbnd) GO TO 220
step = scxbnd

! MOVE SXBD TO THE LEFT SO THAT SBND + TOL(XBND) = XBND.

scxbnd = scxbnd - (reltol*ABS(xbnd)+abstol)/(one + reltol)
220 u = step
IF (ABS(step) < tol .AND. step < zero) u = -tol
IF (ABS(step) < tol .AND. step >= zero) u = tol
itest = 1
RETURN
END SUBROUTINE getptc


FUNCTION dnrm2 ( n, x, incx) RESULT(fn_val)

!  Euclidean norm of the n-vector stored in x() with storage increment incx .
!  if n <= 0 return with result = 0.
!  if n >= 1 then incx must be >= 1

!  c.l.lawson, 1978 jan 08
!  modified to correct failure to update ix, 1/25/92.
!  modified 3/93 to return if incx <= 0.
!  This version by Alan.Miller @ vic.cmis.csiro.au
!  Latest revision - 7 May 2000

!  four phase method using two built-in constants that are
!  hopefully applicable to all machines.
!      cutlo = maximum of  SQRT(u/eps)  over all known machines.
!      cuthi = minimum of  SQRT(v)      over all known machines.
!  where
!      eps = smallest no. such that eps + 1. > 1.
!      u   = smallest positive no.   (underflow limit)
!      v   = largest  no.            (overflow  limit)

!  brief outline of algorithm..

!  phase 1    scans zero components.
!  move to phase 2 when a component is nonzero and <= cutlo
!  move to phase 3 when a component is > cutlo
!  move to phase 4 when a component is >= cuthi/m
!  where m = n for x() real and m = 2*n for complex.

IMPLICIT NONE
INTEGER, INTENT(IN)   :: n, incx
REAL (dp), INTENT(IN) :: x(:)
REAL (dp)             :: fn_val

! Local variables
INTEGER              :: i, ix, j, next
REAL (dp)            :: cuthi, cutlo, hitest, sum, xmax
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp

IF(n <= 0 .OR. incx <= 0) THEN
  fn_val = zero
  RETURN
END IF

! Set machine-dependent constants

cutlo = SQRT( TINY(one) / EPSILON(one) )
cuthi = SQRT( HUGE(one) )

next = 1
sum = zero
i = 1
ix = 1
!                                                 begin main loop
20 SELECT CASE (next)
  CASE (1)
     IF( ABS(x(i)) > cutlo) GO TO 85
     next = 2
     xmax = zero
     GO TO 20

  CASE (2)
!                   phase 1.  sum is zero

     IF( x(i) == zero) GO TO 200
     IF( ABS(x(i)) > cutlo) GO TO 85

!                                prepare for phase 2.   x(i) is very small.
     next = 3
     GO TO 105

  CASE (3)
!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.

     IF( ABS(x(i)) > cutlo ) THEN
!                  prepare for phase 3.

       sum = (sum * xmax) * xmax
       GO TO 85
     END IF

  CASE (4)
     GO TO 110
END SELECT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.

110 IF( ABS(x(i)) <= xmax ) GO TO 115
sum = one + sum * (xmax / x(i))**2
xmax = ABS(x(i))
GO TO 200

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!                   phase 3.  sum is mid-range.  no scaling.

!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)

85 hitest = cuthi / REAL( n, dp )

DO j = ix, n
  IF(ABS(x(i)) >= hitest) GO TO 100
  sum = sum + x(i)**2
  i = i + incx
END DO
fn_val = SQRT( sum )
RETURN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!                                prepare for phase 4.
!                                ABS(x(i)) is very large
100 ix = j
next = 4
sum = (sum / x(i)) / x(i)
!                                Set xmax; large if next = 4, small if next = 3
105 xmax = ABS(x(i))

115 sum = sum + (x(i)/xmax)**2

200 ix = ix + 1
i = i + incx
IF( i <= n ) GO TO 20

!              end of main loop.

!              compute square root and adjust for scaling.

fn_val = xmax * SQRT(sum)

RETURN
END FUNCTION dnrm2

END MODULE trn
