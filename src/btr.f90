MODULE btr
! Global optimization using the Boender-Timmer-Rinnoy Kan algorithm
! The Fortran 77 version was obtained from the web site of Tibor Csendes:
!   www.inf.u-szeged.hu/~csendes/
!
! Modified to interface with FERRE -- C. Allende Prieto, Feb 2011
!

use share, only: dp,ndim,nlambda1,mlsf,nlsf
use fun

IMPLICIT NONE

PRIVATE
PUBLIC :: global

CONTAINS


!   ROUTINE NAME    - GLOBAL
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-09-15  Time: 20:27:48

!-----------------------------------------------------------------------

!   COMPUTER           - IBM PC/SINGLE

!   LATEST REVISION    - OKTOBER 23, 1986

!   PURPOSE            - GLOBAL MINIMUM OF FUNCTION OF N VARIABLES
!                USING A LOCAL SEARCH METHOD

!   USAGE        - CALL GLOBAL (AMIN, AMAX, NPARM, M, N100, NG0, IPR,
!                NSIG, X0, NC, F0)

!   ARGUMENTS   AMIN    - VECTOR OF LENGTH NPARM CONTAINING THE LOWER
!                         BOUNDS OF THE PARAMETERS, SO X(I) IS
!                SEARCHED IN THE INTERVAL (AMIN(I), AMAX(I)).
!                (INPUT)
!        AMAX    - VECTOR OF LENGTH NPARM CONTAINING THE UPPER
!                BOUNDS OF THE PARAMETERS. (INPUT)
!        NPARM   - NUMBER OF PARAMETERS <= 15. (INPUT)
!        M       - NUMBER OF RESIDUAL FUNCTIONS,    WHEN THE
!                OBJECTIVE FUNCTION IS OF THE FORM F1**2+
!                F2**2+...+FM**2, <= 100.    (INPUT)
!        N100    - NUMBER OF SAMPLE POINTS TO BE DRAWN UNIFORMLY
!                IN ONE CYCLE, <= 10000.  THE SUGGESTED VALUE
!                IS 100*NPARM. (INPUT)
!         NG0    - NUMBER OF BEST POINTS SELECTED FROM THE ACTUAL SAMPLE.
!                THE SUGGESTED VALUE IS
!                TWICE THE EXPECTED NUMBER OF LOCAL MINIMA.
!                (INPUT)
!         IPR    - FORTRAN DATA SET REFERENCE NUMBER WHERE THE
!                PRINTED OUTPUT BE SENT. (INPUT)
! ---> changed to verbose level (CAP)
!        NSIG    - CONVERGENCE CRITERION, THE ACCURACY REQUIRED IN THE
!                PARAMETER ESTIMATES.  THIS CONVERGENCE CRITERION IS SATISFIED
!                IF ON TWO SUCCESSIVE ITERATIONS THE PARAMETER
!                ESTIMATES AGREE,
!                COMPONENT BY COMPONENT, TO NSIG DIGITS.
!                THE SUGGESTED VALUE IS 6. (INPUT)
!           X0    - OUTPUT 15 BY 20 MATRIX CONTAINING NC (UP TO 20)
!                   LOCAL MINIMIZERS FOUND.
!           NC    - NUMBER OF DIFFERENT LOCAL MINIMIZERS FOUND.
!                (OUTPUT)
!           F0    - OUTPUT VECTOR OF NC (UP TO 20) OBJECTIVE
!                FUNCTION VALUES, F0(I) BELONGS TO THE
!                PARAMETERS X0(1,I), X0(2,I), ..., X0(NPARM,I).

!   REQUIRED ROUTINES    - URDMN, FUN, LOCAL

!-----------------------------------------------------------------------

SUBROUTINE global(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr,  & 
                  amin, amax, nparm, m, n100, ng0, ipr, nsig, x0, nc, f0)
IMPLICIT NONE


real(dp), intent(in)	:: w(nlambda1)		! weights
real(dp), intent(inout)	:: pf(ndim)			! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)	:: obs(nlambda1)	! vector of observations
real(dp), intent(in)	:: lambda_obs(nlambda1)	! vector of wavelengths
real(dp), intent(in)	:: e_obs(nlambda1)	! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray


REAL(dp), INTENT(IN)         :: amin(15)
REAL(dp), INTENT(IN)         :: amax(15)
INTEGER, INTENT(IN)      :: nparm
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN OUT)  :: n100
INTEGER, INTENT(IN OUT)  :: ng0
INTEGER, INTENT(IN)      :: ipr
INTEGER, INTENT(IN)      :: nsig
REAL(dp), INTENT(IN OUT)     :: x0(15,20)
INTEGER, INTENT(OUT)     :: nc
REAL(dp), INTENT(IN OUT)     :: f0(20)

REAL(dp)     :: x(15,100), x1(15,20), xcl(15,100)
INTEGER  :: ic(100), ic1(20)
REAL(dp)     :: r(100,15), wg(15)
REAL(dp)     :: f(100), f1(20), fcl(100), y(15), mmin(15), mmax(15), b
REAL(dp), PARAMETER  :: zero = 0.0, one = 1.0, two = 2.0, ten = 10.0
INTEGER  :: i, i1, icc, icj, ig, ii, iii, im, in1, inum, inum1, inum2, it, iv, &
            j, jj, l1, maxfn, n, n0, n1, ncp, nfe, nfe1, ng, ng10, nm, nn100, ns
REAL(dp)     :: a, alfa, b1, bb, fc, ff, fm, relcon

IF (nparm <= 0) GO TO 650
IF (nparm > 15) GO TO 640
IF (m > 100) GO TO 650
DO  i = 1, nparm
  mmin(i) = amin(i)
  mmax(i) = amax(i)
  IF (mmin(i) == mmax(i)) GO TO 650
END DO
b1 = one / REAL(nparm)
IF (ng0 < 1) ng0 = 1
IF (ng0 > 20) ng0 = 20
IF (n100 < 20) n100 = 20
IF (n100 > 10000) n100 = 10000
IF (n100 < 100) THEN
  nn100 = n100
  n = 1
ELSE
  nn100 = 100
  n = n100 / 100
  n100 = n * 100
END IF
ng10 = 100
DO  i = 1, ng10
  f(i) = 9.9E10
  ic(i) = 0
END DO
DO  i = 1, nparm
  mmax(i) = (mmax(i)-mmin(i)) / two
  mmin(i) = mmin(i) + mmax(i)
END DO
alfa = .01
nfe = 0
ng = 0
ns = 0
nc = 0
ncp = 1
n0 = 0
n1 = 0
im = 1
ig = 0
fm = 9.9E10
maxfn = 2000 * nparm
relcon = ten ** (-nsig)

!                  SAMPLING
40 n0 = n0 + n100
nm = n0 - 1
ng = ng + ng0
ns = ns + 1
IF (ns*ng0 > 100) GO TO 660
b = (one - alfa**(one/REAL(nm))) ** b1
bb = 0.1 * b
DO  i1 = 1, n
  CALL urdmn(r,1500)
  DO  j = 1, nn100
    DO  i = 1, nparm
      y(i) = two * r(j,i) - one
    END DO
    CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, y, fc, nparm, m, mmin, mmax)
    IF (fc < fm) THEN
      f(im) = fc
      DO  i = 1, nparm
        x(i,im) = y(i)
      END DO
      IF (im <= ng .AND. ic(im) > 0) ig = ig - 1
      ic(im) = 0
      im = 1
      fm = f(1)
      DO  i = 2, ng10
        IF (f(i) >= fm) THEN
          im = i
          fm = f(i)
        END IF
      END DO
    END IF
  END DO
END DO
nfe = nfe + n100
!WRITE (ipr,5000) n100
if (ipr > -10) WRITE (*,5000) n100

!                 SORTING
inum = ng10 - 1
DO  i = 1, inum
  im = i
  fm = f(i)
  inum1 = i + 1
  DO  j = inum1, ng10
    IF (f(j) < fm) THEN
      im = j
      fm = f(j)
    END IF
  END DO
  IF (im > i) THEN
    a = fm
    DO  j = 1, nparm
      y(j) = x(j,im)
    END DO
    IF (i <= ng .AND. im > ng) THEN
      IF (ic(ng) == 0 .AND. ic(im) > 0) ig = ig + 1
      IF (ic(ng) > 0 .AND. ic(im) == 0) ig = ig - 1
    END IF
    icc = ic(im)
    inum1 = im - i
    DO  j = 1, inum1
      inum2 = im - j
      f(inum2+1) = f(inum2)
      ic(inum2+1) = ic(inum2)
      DO  jj = 1, nparm
        x(jj,inum2+1) = x(jj,inum2)
      END DO
    END DO
    f(i) = a
    DO  j = 1, nparm
      x(j,i) = y(j)
    END DO
    ic(i) = icc
  END IF
END DO
IF (nc > 0) THEN

!              CLUSTERING TO    X*
  DO  iii = 1, nc
    i = 1
    in1 = i
    fcl(i) = f0(iii)
    DO  j = 1, nparm
      xcl(j,i) = x0(j,iii)
    END DO
    DO  j = 1, ng
      IF (ic(j) == iii) THEN
        in1 = in1 + 1
        xcl(1:nparm,in1) = x(1:nparm,j)
      END IF
    END DO
    190 DO  j = 1, ng
      IF (ic(j) == 0) THEN
        IF (fcl(i) < f(j)) THEN
          DO  l1 = 1, nparm
            wg(l1) = ABS(xcl(l1,i)-x(l1,j))
          END DO
          a = zero
          DO  l1 = 1, nparm
            IF (wg(l1) > a) a = wg(l1)
          END DO
          IF (a < b) THEN
            !WRITE (ipr,5100) iii
            if (ipr > -10) WRITE (*,5100) iii
            DO  ii = 1, nparm
              wg(ii) = x(ii,j) * mmax(ii) + mmin(ii)
            END DO
            !WRITE (ipr,5200) f(j), (wg(ii),ii = 1,nparm)
            if (ipr > -10) WRITE (*,5200) f(j), (wg(ii),ii = 1,nparm)
            ig = ig + 1
            IF (ig >= ng) GO TO 550
            in1 = in1 + 1
            fcl(in1) = f(j)
            DO  ii = 1, nparm
              xcl(ii,in1) = x(ii,j)
            END DO
            ic(j) = iii
          END IF
        END IF
      END IF
    END DO
    i = i + 1
    IF (i <= in1) GO TO 190
  END DO
  IF (n1 > 0) THEN

!              CLUSTERING TO    X1
    DO  iii = 1, n1
      i = 1
      in1 = i
      fcl(i) = f1(iii)
      DO  j = 1, nparm
        xcl(j,i) = x1(j,iii)
      END DO
      270 DO  j = 1, ng
        IF (ic(j) == 0) THEN
          IF (fcl(i) < f(j)) THEN
            DO  l1 = 1, nparm
              wg(l1) = ABS(xcl(l1,i)-x(l1,j))
            END DO
            a = zero
            DO  l1 = 1, nparm
              IF (wg(l1) > a) a = wg(l1)
            END DO
            IF (a < b) THEN
              !WRITE (ipr,5100) ic1(iii)
              if (ipr > -10) WRITE (*,5100) ic1(iii)
              DO  ii = 1, nparm
                wg(ii) = x(ii,j) * mmax(ii) + mmin(ii)
              END DO
              !WRITE (ipr,5200) f(j), wg(1:nparm)
              if (ipr > -10) WRITE (*,5200) f(j), wg(1:nparm)
              ig = ig + 1
              IF (ig >= ng) GO TO 550
              in1 = in1 + 1
              fcl(in1) = f(j)
              DO  ii = 1, nparm
                xcl(ii,in1) = x(ii,j)
              END DO
              ic(j) = ic1(iii)
            END IF
          END IF
        END IF
      END DO
      i = i + 1
      IF (i <= in1) GO TO 270
    END DO
  END IF
END IF

!              LOCAL SEARCH
it = 0
DO  i1 = 1, ng
  IF (ic(i1) == 0) THEN
    y(1:nparm) = x(1:nparm,i1)
    ff = f(i1)
    CALL local(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, & 
               m, nparm, relcon, maxfn, y, ff, nfe1, r(:,1), mmin, mmax)
    IF (nc > 0) THEN
      DO  iv = 1, nc
        DO  l1 = 1, nparm
          wg(l1) = ABS(x0(l1,iv) - y(l1))
        END DO
        a = zero
        DO  l1 = 1, nparm
          IF (wg(l1) > a) a = wg(l1)
        END DO
        IF (a < bb) GO TO 380
      END DO
      GO TO 430

!                NEW    SEED-POINT
      380 n1 = n1 + 1
      !WRITE (ipr,5300) iv, nfe1
      if (ipr > -10) WRITE (*,5300) iv, nfe1
      DO  ii = 1, nparm
        wg(ii) = x(ii,i1) * mmax(ii) + mmin(ii)
      END DO
      !WRITE (ipr,5200) ff, (wg(ii),ii = 1,nparm)
      if (ipr > -10) WRITE (*,5200) ff, (wg(ii),ii = 1,nparm)
      IF (ff < f0(iv)) THEN
        !WRITE (ipr,5400) iv, f0(iv), ff
        if (ipr > -10) WRITE (*,5400) iv, f0(iv), ff
        DO  ii = 1, nparm
          wg(ii) = y(ii) * mmax(ii) + mmin(ii)
        END DO
        !WRITE (ipr,5200) ff, (wg(ii),ii = 1,nparm)
        if (ipr > -10) WRITE (*,5200) ff, (wg(ii),ii = 1,nparm)
        f0(iv) = ff
        DO  ii = 1, nparm
          x0(ii,iv) = y(ii)
        END DO
      END IF
      IF (n1 > 20) GO TO 670
      DO  ii = 1, nparm
        x1(ii,n1) = x(ii,i1)
        xcl(ii,1) = x(ii,i1)
      END DO
      f1(n1) = f(i1)
      fcl(1) = f(i1)
      ic1(n1) = iv
      icj = iv
      GO TO 460
    END IF

!              NEW LOCAL MINIMUM
    430 nc = nc + 1
    ncp = ncp + 1
    !WRITE (ipr,5500) nc, ff, nfe1
    if (ipr > -10) WRITE (*,5500) nc, ff, nfe1
    DO  ii = 1, nparm
      wg(ii) = y(ii) * mmax(ii) + mmin(ii)
    END DO
    !WRITE (ipr,5200) ff, (wg(ii),ii = 1,nparm)
    if (ipr > -10) WRITE (*,5200) ff, (wg(ii),ii = 1,nparm)
    DO  ii = 1, nparm
      x0(ii,nc) = y(ii)
      xcl(ii,1) = y(ii)
    END DO
    fcl(1) = ff
    f0(nc) = ff
    IF (nc >= 20) GO TO 680
    it = 1
    icj = nc

!             CLUSTERING TO THE NEW POINT
    460 nfe = nfe + nfe1
    ic(i1) = icj
    ig = ig + 1
    IF (ig >= ng) EXIT
    i = 1
    in1 = i
    470 DO  j = 1, ng
      IF (ic(j) == 0) THEN
        IF (fcl(i) < f(j)) THEN
          DO  l1 = 1, nparm
            wg(l1) = ABS(xcl(l1,i)-x(l1,j))
          END DO
          a = zero
          DO  l1 = 1, nparm
            IF (wg(l1) > a) a = wg(l1)
          END DO
          IF (a < b) THEN
            in1 = in1 + 1
            DO  ii = 1, nparm
              xcl(ii,in1) = x(ii,j)
            END DO
            fcl(in1) = f(j)
            ic(j) = icj
            !WRITE (ipr,5100) icj
            if (ipr > -10) WRITE (*,5100) icj
            DO  ii = 1, nparm
              wg(ii) = x(ii,j) * mmax(ii) + mmin(ii)
            END DO
            !WRITE (ipr,5200) f(j), (wg(ii),ii = 1,nparm)
            if (ipr > -10) WRITE (*,5200) f(j), (wg(ii),ii = 1,nparm)
            ig = ig + 1
            IF (ig >= ng) EXIT
          END IF
        END IF
      END IF
    END DO
    i = i + 1
    IF (i < in1) GO TO 470
  END IF
END DO
IF (it /= 0) GO TO 40

!               PRINT RESULTS
!550 WRITE (ipr,5600)
550 WRITE (*,5600)
IF (nc > 1) THEN
  inum = nc - 1
  DO  i = 1, inum
    im = i
    fm = f0(i)
    inum1 = i + 1
    DO  j = inum1, nc
      IF (f0(j) < fm) THEN
        im = j
        fm = f0(j)
      END IF
    END DO
    IF (im > i) THEN
      a = fm
      DO  j = 1, nparm
        y(j) = x0(j,im)
      END DO
      inum1 = im - i
      DO  j = 1, inum1
        inum2 = im - j
        f0(inum2+1) = f0(inum2)
        DO  jj = 1, nparm
          x0(jj,inum2+1) = x0(jj,inum2)
        END DO
      END DO
      f0(i) = a
      DO  j = 1, nparm
        x0(j,i) = y(j)
      END DO
    END IF
  END DO
END IF

IF (nc > 0) THEN
  DO  i = 1, nc
    DO  ii = 1, nparm
      x0(ii,i) = x0(ii,i) * mmax(ii) + mmin(ii)
    END DO
    !WRITE (ipr,5200) f0(i), (x0(ii,i),ii = 1,nparm)
    WRITE (*,5200) f0(i), (x0(ii,i),ii = 1,nparm)
  END DO
END IF
!WRITE (ipr,5700) nfe
WRITE (*,5700) nfe
RETURN
!640 WRITE (ipr,5800)
640 WRITE (*,5800)
STOP
!650 WRITE (ipr,5900)
650 WRITE (*,5900)
STOP
!660 WRITE (ipr,6000)
660 WRITE (*,6000)
GO TO 550
!670 WRITE (ipr,6100)
670 WRITE (*,6100)
GO TO 550
!680 WRITE (ipr,6200)
680 WRITE (*,6200)
GO TO 550

!STOP

5000 FORMAT (/' ', i5 ,' FUNCTION EVALUATIONS USED FOR SAMPLING')
5100 FORMAT (' SAMPLE POINT ADDED TO THE CLUSTER NO. ',i2)
5200 FORMAT (' ', g15.8, 3(/t5, 5(g15.8, ' ')))
5300 FORMAT (' NEW SEED POINT ADDED TO THE CLUSTER NO. ',i2,', NFEV=', i5)
5400 FORMAT (' *** IMPROVEMENT ON THE LOCAL MINIMUM NO. ',i2,':',g15.8,  &
    ' FOR ',g15.8)
5500 FORMAT (' *** THE LOCAL MINIMUM NO. ',i2,': ',g15.8,', NFEV=',i5)
5600 FORMAT (/////,' LOCAL MINIMA FOUND:'//)
5700 FORMAT (///,' NORMAL TERMINATION AFTER ',i5,' FUNCTION ',  &
    'EVALUATIONS',///)
5800 FORMAT (' ***   TOO MANY PARAMETERS',//,' ABNORMAL TERMINATION')
5900 FORMAT (' ***   DATA ERROR')
6000 FORMAT (' ***   TOO MANY SAMPLING')
6100 FORMAT (' ***   TOO MANY NEW SEED POINTS')
6200 FORMAT (' ***   TOO MANY CLUSTERS')
END SUBROUTINE global



SUBROUTINE globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, r, f, nparm, m, MIN, MAX)
IMPLICIT NONE

real(dp), intent(in)	:: w(nlambda1)		! weights
real(dp), intent(inout)	:: pf(ndim)			! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)	:: obs(nlambda1)	! vector of observations
real(dp), intent(in)	:: lambda_obs(nlambda1)	! vector of wavelengths
real(dp), intent(in)	:: e_obs(nlambda1)	! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray

REAL(dp), INTENT(IN)      :: r(:)
REAL(dp), INTENT(IN OUT)  :: f
INTEGER, INTENT(IN)   :: nparm
INTEGER, INTENT(IN)   :: m
REAL(dp), INTENT(IN)      :: MIN(15)
REAL(dp), INTENT(IN)      :: MAX(15)


REAL(dp)     :: x(15)
INTEGER  :: i

DO  i = 1, nparm
  x(i) = MAX(i) * r(i) + MIN(i)
END DO
CALL objfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, x, f)
RETURN

END SUBROUTINE globalfun



!   ROUTINE NAME    - LOCAL
 
!-----------------------------------------------------------------------

!   COMPUTER            - IBM PC/SINGLE

!   LATEST REVISION    - JULY 31, 1986

!   PURPOSE        - MINIMUM OF A FUNCTION OF N VARIABLES USING
!                  A QUASI-NEWTON METHOD

!   USAGE        - CALL LOCAL (M,N,EPS,MAXFN,X,F,NFEV,W,MIN,MAX)

!   ARGUMENTS
!         M    - THE NUMBER OF RESIDUAL FUNCTIONS (INPUT)
!                NOT USED IN THIS ROUTINE.
!         N    - THE NUMBER OF PARAMETERS (I.E., THE LENGTH
!                OF X) (INPUT)
!         EPS    - CONVERGENCE CRITERION. (INPUT). THE ACCURACY
!                REQUIRED IN THE PARAMETER ESTIMATES
!                THIS CONVERGENCE CONDITION IS SATISFIED IF
!                ON TWO SUCCESSIVE ITERATIONS, THE PARAMETER
!                ESTIMATES (I.E.,X(I), I=1,...,N) DIFFERS,
!                COMPONENT BY COMPONENT, BY AT MOST EPS.
!         MAXFN    - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E.,
!                CALLS TO SUBROUTINE FUN) ALLOWED. (INPUT)
!         X    - VECTOR OF LENGTH N CONTAINING PARAMETER VALUES.
!              ON INPUT, X MUST CONTAIN THE INITIAL
!                PARAMETER ESTIMATES.
!              ON OUTPUT, X CONTAINS    THE FINAL PARAMETER
!                ESTIMATES AS DETERMINED BY LOCAL.
!         F    - A SCALAR CONTAINING THE VALUE    OF THE FUNCTION
!                AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)
!         NFEV    - THE NUMBER OF FUNCTION EVALUATIONS (OUTPUT)
!         W    - A VECTOR OF LENGTH 3*N USED AS WORKING SPACE.
!                MMIN    - A VECTOR OF LENGTH N CONTAINING THE LOWER
!                           BOUNDS OF THE PARAMETERS, SO X(I) IS
!                           SEARCHED IN THE INTERVAL (MIN(I),MAX(I)).
!                           (INPUT)
!                MMAX    - A VECTOR OF LENGTH N CONTAINING THE UPPER
!                           BOUNDS OF THE PARAMETERS. (INPUT)

!   REQUIRED ROUTINES    - UPDATE, FUN

!         FUN      - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
!                THE FUNCTION F FOR GIVEN PARAMETER VALUES
!                X(1),X(2),...,X(N).
!                THE CALLING SEQUENCE HAS THE FOLLOWING FORM
!                CALL FUN(X, F, N, M, MIN, MAX)
!                WHERE X IS A VEKTOR OF LENGTH N.
!                FUN MUST APPEAR IN AN EXTERNAL STATEMENT
!                IN THE CALLING PROGRAM.  FUN MUST NOT
!                ALTER THE VALUES OF    X(I),I=1,...,N OR N.

!-----------------------------------------------------------------------

SUBROUTINE local(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, & 
                 m, n, eps, maxfn, x, f, nfev, wg, mmin, mmax)
!                   SPECIFICATIONS FOR ARGUMENTS

IMPLICIT NONE

real(dp), intent(in)	:: w(nlambda1)		! weights
real(dp), intent(inout)	:: pf(ndim)			! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)	:: obs(nlambda1)	! vector of observations
real(dp), intent(in)	:: lambda_obs(nlambda1)	! vector of wavelengths
real(dp), intent(in)	:: e_obs(nlambda1)	! vector of uncertainties
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray

INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN)      :: n
REAL(dp), INTENT(IN)         :: eps
INTEGER, INTENT(IN OUT)  :: maxfn
REAL(dp), INTENT(IN OUT)     :: x(n)
REAL(dp), INTENT(IN OUT)     :: f
INTEGER, INTENT(OUT)     :: nfev
REAL(dp), INTENT(OUT)        :: wg(:)
REAL(dp), INTENT(IN OUT)     :: mmin(:)
REAL(dp), INTENT(IN OUT)     :: mmax(:)

!                   SPECIFICATIONS FOR LOCAL VARIABLES
INTEGER :: ig, igg, is, idiff, ir, ij, i, iopt, j, nm1, jj, jp1, l, kj, k,  &
           link, itn, ii, im1, jnt, np1, jb, nj, ier
REAL(dp)    :: hh, hjj, v, df, relx, gs0, diff, aeps, alpha, ff, tot, f1, f2,  &
           z, gys, dgs, sig, zz, hhh, ghh, g(15), h(120)
REAL(dp), PARAMETER  :: reps = 1.1921E-07, zero = 0.0, one = 1.0, half = 0.5, &
                    seven = 7.0, five = 5.0, twelve = 12.0, p1 = 0.1

!                   INITIALIZATION
!                   FIRST EXECUTABLE STATEMENT
iopt = 0
!         IOPT    - OPTIONS SELECTOR. (INPUT)
!              IOPT = 0 CAUSES LOCAL TO INITIALIZE THE
!                HESSIAN MATRIX H TO THE IDENTITY MATRIX.
!              IOPT = 1 INDICATES THAT H HAS BEEN INITIALIZED
!                BY THE USER TO A POSITIVE DEFINITE MATRIX.
!              IOPT = 2 CAUSES LOCAL TO COMPUTE THE DIAGONAL
!                VALUES OF THE HESSIAN MATRIX AND SET H TO
!                A DIAGONAL MATRIX CONTAINING THESE VALUES.
!              IOPT = 3 CAUSES LOCAL TO COMPUTE AN ESTIMATE
!                OF THE HESSIAN IN H.
ier = 0
hh = SQRT(reps)
ig = n
igg = n + n
is = igg
idiff = 1
ir = n
wg(1) = -one
wg(2) = zero
wg(3) = zero

!                   EVALUATE FUNCTION AT STARTING POINT
g(1:n) = x(1:n)
CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, g, f, n, m, mmin, mmax)
nfev = 1
IF (iopt /= 1) THEN
!                   SET OFF-DIAGONAL ELEMENTS OF H TO 0.0
  IF (n /= 1) THEN
    ij = 2
    DO  i = 2, n
      DO  j = 2, i
        h(ij) = zero
        ij = ij + 1
      END DO
      ij = ij + 1
    END DO
    IF (iopt == 0) THEN
!                   SET DIAGONAL ELEMENTS OF H TO ONE
      ij = 0
      DO  i = 1, n
        ij = ij + i
        h(ij) = one
      END DO
      GO TO 110
    END IF
  END IF
!                   GET DIAGONAL ELEMENTS OF HESSIAN
  im1 = 1
  nm1 = 1
  np1 = n + 1
  DO  i = 2, np1
    hhh = hh * ABS(x(im1))
    g(im1) = x(im1) + hhh
    CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, g, f2, n, m, mmin, mmax)
    g(im1) = g(im1) + hhh
    CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, g, ff, n, m, mmin, mmax)
    h(nm1) = (ff-f2+f-f2) / (hhh*hhh)
    g(im1) = x(im1)
    im1 = i
    nm1 = i + nm1
  END DO
  nfev = nfev + n + n
  IF (iopt == 3 .AND. n /= 1) THEN
!                   GET THE REST OF THE HESSIAN
    jj = 1
    ii = 2
    DO  i = 2, n
      ghh = hh * ABS(x(i))
      g(i) = x(i) + ghh
      CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, g, f2, n, m, mmin, mmax)
      DO  j = 1, jj
        hhh = hh * ABS(x(j))
        g(j) = x(j) + hhh
        CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, g, ff, n, m, mmin, mmax)
        g(i) = x(i)
        CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, g, f1, n, m, mmin, mmax)
!        H(II) = (FF-F1-F2+F)*SQREPS
        h(ii) = (ff-f1-f2+f) / (hhh*ghh)
        ii = ii + 1
        g(j) = x(j)
      END DO
      jj = jj + 1
      ii = ii + 1
    END DO
    nfev = nfev + ((n*n-n)/2)
  END IF
END IF
!                   FACTOR H TO L*D*L-TRANSPOSE
ir = n
IF (n <= 1) THEN
  IF (h(1) > zero) GO TO 110
  h(1) = zero
  ir = 0
ELSE
  nm1 = n - 1
  jj = 0
  DO  j = 1, n
    jp1 = j + 1
    jj = jj + j
    hjj = h(jj)
    IF (hjj <= zero) THEN
      h(jj) = zero
      ir = ir - 1
    ELSE
      IF (j /= n) THEN
        ij = jj
        l = 0
        DO  i = jp1, n
          l = l + 1
          ij = ij + i - 1
          v = h(ij) / hjj
          kj = ij
          DO  k = i, n
            h(kj+l) = h(kj+l) - h(kj) * v
            kj = kj + k
          END DO
          h(ij) = v
        END DO
      END IF
    END IF
  END DO
END IF
IF (ir /= n) THEN
  ier = 129
  GO TO 440
END IF
110 itn = 0
df = -one

!                   EVALUATE GRADIENT W(IG+I),I=1,...,N
120 link = 1
GO TO 410

!                   BEGIN ITERATION LOOP
130 IF (nfev >= maxfn) THEN
    ier=200 ! Added to avoid an infinite loop (C. Allende Prieto, April 4, 2017)
    write(*,*)'Maximum number of function evaluations exceeded. Giving up!'
    GO TO 350
ENDIF
itn = itn + 1
DO  i = 1, n
  wg(i) = -wg(ig+i)
END DO

!                   DETERMINE SEARCH DIRECTION W
!                     BY    SOLVING    H*W = -G WHERE
!                     H = L*D*L-TRANSPOSE
IF (ir >= n) THEN
!                   N .EQ. 1
  g(1) = wg(1)
  IF (n <= 1) THEN
    wg(1) = wg(1) / h(1)
  ELSE
!                   N .GT. 1
    ii = 1
!                   SOLVE L*W = -G
    DO  i = 2, n
      ij = ii
      ii = ii + i
      v = wg(i)
      im1 = i - 1
      DO  j = 1, im1
        ij = ij + 1
        v = v - h(ij) * wg(j)
      END DO
      g(i) = v
      wg(i) = v
    END DO
!                   SOLVE (D*LT)*Z = W WHERE LT = L-TRANSPOSE
    wg(n) = wg(n) / h(ii)
    jj = ii
    nm1 = n - 1
    DO  nj = 1, nm1
!                   J = N-1,N-2,...,1
      j = n - nj
      jp1 = j + 1
      jj = jj - jp1
      v = wg(j) / h(jj)
      ij = jj
      DO  i = jp1, n
        ij = ij + i - 1
        v = v - h(ij) * wg(i)
      END DO
      wg(j) = v
    END DO
  END IF
END IF

!                   DETERMINE STEP LENGTH ALPHA
relx = zero
gs0 = zero
DO  i = 1, n
  wg(is+i) = wg(i)
  diff = ABS(wg(i)) / ABS(x(i))
  relx = MAX(relx, diff)
  gs0 = gs0 + wg(ig+i) * wg(i)
END DO
IF (relx == zero) GO TO 360
aeps = eps / relx
ier = 130
IF (gs0 >= zero) GO TO 360
IF (df == zero) GO TO 360
ier = 0
alpha = (-df-df) / gs0
IF (alpha <= zero) alpha = one
alpha = MIN(alpha, one)
IF (idiff == 2) alpha = MAX(p1,alpha)
ff = f
tot = zero
jnt = 0

!                   SEARCH ALONG    X + ALPHA*W
200 IF (nfev >= maxfn) GO TO 350
DO  i = 1, n
  wg(i) = x(i) + alpha * wg(is+i)
END DO
CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, wg, f1, n, m, mmin, mmax)
nfev = nfev + 1
IF (f1 < f) THEN
  f2 = f
  tot = tot + alpha
  220 ier = 0
  f = f1
  DO  i = 1, n
    x(i) = wg(i)
  END DO
  IF (jnt-1 < 0) THEN
    GO TO   240
  ELSE IF (jnt-1 == 0) THEN
    GO TO   280
  ELSE
    GO TO   290
  END IF
  240 IF (nfev >= maxfn) GO TO 350
  DO  i = 1, n
    wg(i) = x(i) + alpha * wg(is+i)
  END DO
  CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, wg, f1, n, m, mmin, mmax)
  nfev = nfev + 1
  IF (f1 >= f) GO TO 290
  IF (f1+f2 >= f+f .AND. seven*f1+five*f2 > twelve*f) jnt = 2
  tot = tot + alpha
  alpha = alpha + alpha
  GO TO 220
END IF
IF (f == ff .AND. idiff == 2 .AND. relx > eps) ier = 130
IF (alpha < aeps) GO TO 360
IF (nfev >= maxfn) GO TO 350
alpha = half * alpha
DO  i = 1, n
  wg(i) = x(i) + alpha * wg(is+i)
END DO
CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, wg, f2, n, m, mmin, mmax)
nfev = nfev + 1
IF (f2 < f) THEN
  tot = tot + alpha
  ier = 0
  f = f2
  DO  i = 1, n
    x(i) = wg(i)
  END DO
ELSE
  z = p1
  IF (f1+f > f2+f2) z = one + half * (f-f1) / (f+f1-f2-f2)
  z = MAX(p1,z)
  alpha = z * alpha
  jnt = 1
  GO TO 200
END IF
280 IF (tot < aeps) GO TO 360
290 alpha = tot

!                   SAVE OLD GRADIENT
DO  i = 1, n
  wg(i) = wg(ig+i)
END DO
!                   EVALUATE GRADIENT W(IG+I), I=1,...,N
link = 2
GO TO 410
310 IF (nfev < maxfn) THEN
  gys = zero
  DO  i = 1, n
    gys = gys + wg(ig+i) * wg(is+i)
    wg(igg+i) = wg(i)
  END DO
  df = ff - f
  dgs = gys - gs0
  IF (dgs <= zero) GO TO 130
  IF (dgs+alpha*gs0 <= zero) THEN
!                   UPDATE HESSIAN H USING COMPLEMENTARY DFP FORMULA
    sig = one / gs0
    ir = -ir
    CALL update(h, n, wg, sig, g, ir, 0, zero)
    DO  i = 1, n
      g(i) = wg(ig+i) - wg(igg+i)
    END DO
    sig = one / (alpha*dgs)
    ir = -ir
    CALL update(h, n, g, sig, wg, ir, 0, zero)
    GO TO 130
  END IF

!                   UPDATE HESSIAN USING DFP FORMULA
  zz = alpha / (dgs-alpha*gs0)
  sig = -zz
  CALL update(h, n, wg, sig, g, ir, 0, reps)
  z = dgs * zz - one
  DO  i = 1, n
    g(i) = wg(ig+i) + z * wg(igg+i)
  END DO
  sig = one / (zz*dgs*dgs)
  CALL update(h, n, g, sig, wg, ir, 0, zero)
  GO TO 130
END IF

!                   MAXFN FUNCTION EVALUATIONS
350 GO TO 370
360 IF (idiff /= 2) THEN
!                   CHANGE TO CENTRAL DIFFERENCES
  idiff = 2
  GO TO 120
END IF
370 IF (relx > eps .AND. ier == 0) GO TO 120

!                   COMPUTE H = L*D*L-TRANSPOSE AND OUTPUT
IF (n == 1) GO TO 440
np1 = n + 1
nm1 = n - 1
jj = (n*(np1)) / 2
DO  jb = 1, nm1
  jp1 = np1 - jb
  jj = jj - jp1
  hjj = h(jj)
  ij = jj
  l = 0
  DO  i = jp1, n
    l = l + 1
    ij = ij + i - 1
    v = h(ij) * hjj
    kj = ij
    DO  k = i, n
      h(kj+l) = h(kj+l) + h(kj) * v
      kj = kj + k
    END DO
    h(ij) = v
  END DO
  hjj = h(jj)
END DO
GO TO 440

!                    EVALUATE GRADIENT
410 IF (idiff /= 2) THEN
!                   FORWARD DIFFERENCES
!                     GRADIENT =    W(IG+I), I=1,...,N
  DO  i = 1, n
    z = hh * ABS(x(i))
    zz = x(i)
    x(i) = zz + z
    CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, x, f1, n, m, mmin, mmax)
    wg(ig+i) = (f1-f) / z
    x(i) = zz
  END DO
  nfev = nfev + n
  SELECT CASE ( link )
    CASE (    1)
      GO TO 130
    CASE (    2)
      GO TO 310
  END SELECT
END IF
!                   CENTRAL DIFFERENCES
!                     GRADIENT =    W(IG+I), I=1,...,N
DO  i = 1, n
  z = hh * ABS(x(i))
  zz = x(i)
  x(i) = zz + z
  CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, x, f1, n, m, mmin, mmax)
  x(i) = zz - z
  CALL globalfun(w, pf, pf0, obs, lambda_obs, e_obs, mobs, lsfarr, x, f2, n, m, mmin, mmax)
  wg(ig+i) = (f1-f2) / (z+z)
  x(i) = zz
END DO
nfev = nfev + n + n
SELECT CASE ( link )
  CASE (    1)
    GO TO 130
  CASE (    2)
    GO TO 310
END SELECT
!                   RETURN
440 RETURN
END SUBROUTINE local



!   ROUTINE NAME    - UPDATE

!-----------------------------------------------------------------------

!   COMPUTER        - IBM PC/SINGLE

!   LATEST REVISION     - JULY 31, 1986
!              (CHANGES IN COMMENTS)

!   PURPOSE        - NUCLEUS CALLED ONLY BY ROUTINE LOCAL

!   PRECISION/HARDWARE    - SINGLE AND DOUBLE/H32
!            - DOUBLE/H36,H48,H60

!   REQD. ROUTINES    - NONE REQUIRED

!-----------------------------------------------------------------------

SUBROUTINE update(a, n, z, sig, wg, ir, mk, eps)
!                   SPECIFICATIONS FOR ARGUMENTS
IMPLICIT NONE

REAL(dp), INTENT(OUT)     :: a(:)
INTEGER, INTENT(IN)   :: n
REAL(dp), INTENT(IN OUT)  :: z(n)
REAL(dp), INTENT(IN)      :: sig
REAL(dp), INTENT(IN OUT)  :: wg(n)
INTEGER, INTENT(OUT)  :: ir
INTEGER, INTENT(IN)   :: mk
REAL(dp), INTENT(IN)      :: eps


!                   SPECIFICATIONS FOR LOCAL VARIABLES
INTEGER :: j, jj, ij, jp1, i, ii, mm
REAL(dp)    :: ti, v, tim, al, r, b, gm, y
REAL(dp), PARAMETER  :: zero = 0.0, one = 1.0, four = 4.0

!                   UPDATE FACTORS GIVEN    IN A
!                     SIG*Z*Z-TRANSPOSE IS ADDED
!                   FIRST EXECUTABLE STATEMENT
IF (n <= 1) THEN
!                   N .EQ. 1
  a(1) = a(1) + sig * z(1) * z(1)
  ir = 1
  IF (a(1) > zero) GO TO 150
  a(1) = zero
  ir = 0
ELSE
!                   N .GT. 1
  IF (sig <= zero) THEN
    IF (sig == zero .OR. ir == 0) GO TO 150
    ti = one / sig
    jj = 0
    IF (mk /= 0) THEN
!                   L*W = Z ON INPUT
      DO  j = 1, n
        jj = jj + j
        IF (a(jj) /= zero) ti = ti + (wg(j)*wg(j)) / a(jj)
      END DO
    ELSE
!                   SOLVE L*W = Z
      wg(1:n) = z(1:n)
      DO  j = 1, n
        jj = jj + j
        v = wg(j)
        IF (a(jj) <= zero) THEN
          wg(j) = zero
        ELSE
          ti = ti + (v*v) / a(jj)
          IF (j /= n) THEN
            ij = jj
            jp1 = j + 1
            DO  i = jp1, n
              ij = ij + i - 1
              wg(i) = wg(i) - v * a(ij)
            END DO
          END IF
        END IF
      END DO
    END IF
!                    SET    TI, TIM    AND W
    IF (ir > 0) THEN
      IF (ti > zero) GO TO 50
      IF (mk-1 > 0) THEN
        GO TO 60
      ELSE
        GO TO 80
      END IF
    END IF
    ti = zero
    ir = -ir - 1
    GO TO 60
    50 ti = eps / sig
    IF (eps == zero) ir = ir - 1
    60 tim = ti
    ii = jj
    i = n
    DO  j = 1, n
      IF (a(ii) /= zero) tim = ti - (wg(i)*wg(i)) / a(ii)
      wg(i) = ti
      ti = tim
      ii = ii - i
      i = i - 1
    END DO
    mm = 1
    GO TO 90
  END IF
  80 mm = 0
  tim = one / sig
  90 jj = 0
!                   UPDATE A
  DO  j = 1, n
    jj = jj + j
    ij = jj
    jp1 = j + 1
!                   UPDATE A(J,J)
    v = z(j)
    IF (a(jj) <= zero) THEN
!                   A(J,J) .EQ. ZERO
      IF (ir <= 0 .AND. sig >= zero .AND. v /= zero) THEN
        ir = 1 - ir
        a(jj) = (v*v) / tim
        IF (j == n) GO TO 150
        DO  i = jp1, n
          ij = ij + i - 1
          a(ij) = z(i) / v
        END DO
        GO TO 150
      END IF
      ti = tim
    ELSE
!                   A(J,J) .GT. ZERO
      al = v / a(jj)
      ti = wg(j)
      IF (mm == 0) ti = tim + v * al
      r = ti / tim
      a(jj) = r * a(jj)
      IF (r == zero) EXIT
      IF (j == n) EXIT
!                   UPDATE REMAINDER OF COLUMN J
      b = al / ti
      IF (r <= four) THEN
        DO  i = jp1, n
          ij = ij + i - 1
          z(i) = z(i) - v * a(ij)
          a(ij) = a(ij) + b * z(i)
        END DO
      ELSE
        gm = tim / ti
        DO  i = jp1, n
          ij = ij + i - 1
          y = a(ij)
          a(ij) = b * z(i) + y * gm
          z(i) = z(i) - v * y
        END DO
      END IF
      tim = ti
    END IF
  END DO
  IF (ir < 0) ir = -ir
END IF

150 RETURN
END SUBROUTINE update



SUBROUTINE urdmn(x, n)
 
! Generates n random numbers
! This replaces the original routine

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL(dp), INTENT(OUT)    :: x(n)

CALL RANDOM_NUMBER(x)
RETURN
END SUBROUTINE urdmn

END MODULE btr
