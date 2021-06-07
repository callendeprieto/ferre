module nelder_mead

use share, only: dp,ndim,nlambda1,mlsf,nlsf
use fun

IMPLICIT NONE

PRIVATE
PUBLIC :: minim

CONTAINS

!*****************************************************************************
SUBROUTINE minim(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p, & 
                 step, nop, func, maxfn, &
		 iprint, stopcr, nloop, iquad, simp, var, ifault)

!     A PROGRAM FOR FUNCTION MINIMIZATION USING THE SIMPLEX METHOD.

!     FOR DETAILS, SEE NELDER & MEAD, THE COMPUTER JOURNAL, JANUARY 1965

!     PROGRAMMED BY D.E.SHAW,
!     CSIRO, DIVISION OF MATHEMATICS & STATISTICS
!     P.O. BOX 218, LINDFIELD, N.S.W. 2070

!     WITH AMENDMENTS BY R.W.M.WEDDERBURN
!     ROTHAMSTED EXPERIMENTAL STATION
!     HARPENDEN, HERTFORDSHIRE, ENGLAND

!     Further amended by Alan Miller
!     CSIRO Division of Mathematical & Information Sciences
!     Private Bag 10, CLAYTON, VIC. 3169

!     Fortran 90 conversion by Alan Miller, June 1995
!     Alan.Miller @ vic.cmis.csiro.au
!     Latest revision - 5 December 1999

!	  Adapted for FERRE, Carlos Allende Prieto, 2002
!

!     ARGUMENTS:-
!     P()     = INPUT, STARTING VALUES OF PARAMETERS
!               OUTPUT, FINAL VALUES OF PARAMETERS
!     STEP()  = INPUT, INITIAL STEP SIZES
!     NOP     = INPUT, NO. OF PARAMETERS, INCL. ANY TO BE HELD FIXED
!     FUNC    = OUTPUT, THE FUNCTION VALUE CORRESPONDING TO THE FINAL
!                 PARAMETER VALUES.
!     maxfn     = INPUT, THE MAXIMUM NO. OF FUNCTION EVALUATIONS ALLOWED.
!               Say, 20 times the number of parameters, NOP.
!     IPRINT  = INPUT, PRINT CONTROL PARAMETER
!                 < 0 NO PRINTING
!                 = 0 PRINTING OF PARAMETER VALUES AND THE FUNCTION
!                     VALUE AFTER INITIAL EVIDENCE OF CONVERGENCE.
!                 > 0 AS FOR IPRINT = 0 PLUS PROGRESS REPORTS AFTER EVERY
!                     IPRINT EVALUATIONS, PLUS PRINTING FOR THE INITIAL SIMPLEX.
!     STOPCR  = INPUT, STOPPING CRITERION.
!               The criterion is applied to the standard deviation of
!               the values of FUNC at the points of the simplex.
!     NLOOP   = INPUT, THE STOPPING RULE IS APPLIED AFTER EVERY NLOOP
!               FUNCTION EVALUATIONS.   Normally NLOOP should be slightly
!               greater than NOP, say NLOOP = 2*NOP.
!     IQUAD   = INPUT, = 1 IF FITTING OF A QUADRATIC SURFACE IS REQUIRED
!                      = 0 IF NOT
!               N.B. The fitting of a quadratic surface is strongly
!               recommended, provided that the fitted function is
!               continuous in the vicinity of the minimum.   It is often
!               a good indicator of whether a premature termination of
!               the search has occurred.
!     SIMP    = INPUT, CRITERION FOR EXPANDING THE SIMPLEX TO OVERCOME
!               ROUNDING ERRORS BEFORE FITTING THE QUADRATIC SURFACE.
!               The simplex is expanded so that the function values at
!               the points of the simplex exceed those at the supposed
!               minimum by at least an amount SIMP.
!     VAR()   = OUTPUT, CONTAINS THE DIAGONAL ELEMENTS OF THE INVERSE OF
!               THE INFORMATION MATRIX.
!     FUNCTN  = INPUT, NAME OF THE USER'S SUBROUTINE - ARGUMENTS (P,FUNC)
!               WHICH RETURNS THE FUNCTION VALUE FOR A GIVEN SET OF
!               PARAMETER VALUES IN ARRAY P.
!****     FUNCTN MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM.
!     IFAULT  = OUTPUT, = 0 FOR SUCCESSFUL TERMINATION
!                 = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
!                 = 2 IF INFORMATION MATRIX IS NOT +VE SEMI-DEFINITE
!                 = 3 IF NOP < 1
!                 = 4 IF NLOOP < 1

!     N.B. P, STEP AND VAR (IF IQUAD = 1) MUST HAVE DIMENSION AT LEAST NOP
!          IN THE CALLING PROGRAM.

!*****************************************************************************

INTEGER, INTENT(IN)        :: nop, maxfn, iprint, nloop, iquad
INTEGER, INTENT(OUT)       :: ifault
REAL (dp), INTENT(IN)      :: stopcr, simp
real(dp), intent(in)	   :: w(nlambda1)			! weights
real(dp), intent(in)    :: chiscale		!chi**2/sum(w_i(f_i-obs_i))^2)
real(dp), intent(inout)	   :: pf(ndim)			! vector of fixed parameters
real(dp), intent(in)       :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)	   :: obs(nlambda1)		! vector of observations
real(dp), intent(in)	   :: lambda_obs(nlambda1)	! vector of wavelengths
real(dp), intent(in)	   :: e_obs(nlambda1)	! vector of uncertainties
real(dp), intent(in)       :: mobs             ! mean or median of obs array
real(dp), intent(in)       :: lsfarr(mlsf,nlsf)    ! lsfarray
REAL (dp), INTENT(IN OUT)  :: p(:), step(:)
REAL (dp), INTENT(OUT)     :: var(:), func

!     Local variables

REAL (dp)   :: g(nop+1,nop), h(nop+1), pbar(nop), pstar(nop), pstst(nop), &
               aval(nop), pmin(nop), temp(nop), bmat(nop*(nop+1)/2),  &
               vc(nop*(nop+1)/2), ymin, rmax, hstst, a0, hmin, test, hmean, &
               hstd, hstar, hmax, savemn

REAL (dp), PARAMETER :: zero = 0._dp, half = 0.5_dp, one = 1._dp, two = 2._dp
INTEGER     :: i, i1, i2, iflag, ii, ij, imax, imin, irank, irow, j, j1, jj, &
               k, l, loop, nap, neval, nmore, np1, nullty

!     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
!     C = EXPANSION COEFFICIENT.

REAL (dp), PARAMETER :: a = 1._dp, b = 0.5_dp, c = 2._dp

!     SET LOUT = LOGICAL UNIT NO. FOR OUTPUT

INTEGER, PARAMETER :: lout = 6

!     IF PROGRESS REPORTS HAVE BEEN REQUESTED, PRINT HEADING

IF (iprint > 0) WRITE (lout,5000) iprint

!     CHECK INPUT ARGUMENTS

ifault = 0
IF (nop <= 0) ifault = 3
IF (nloop <= 0) ifault = 4
IF (ifault /= 0) RETURN

!     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP /= 0

nap = COUNT(step /= zero)
neval = 0
loop = 0
iflag = 0

!     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN

IF (nap <= 0) THEN
  CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,func)
  RETURN
END IF

!     SET UP THE INITIAL SIMPLEX

20 g(1,:) = p
irow = 2
DO i = 1, nop
  IF (step(i) /= zero) THEN
    g(irow,:) = p
    g(irow,i) = p(i) + step(i)
    irow = irow + 1
  END IF
END DO

np1 = nap + 1
DO i = 1, np1
  p = g(i,:)
  CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,h(i))
  neval = neval + 1
  IF (iprint > 0) THEN
    WRITE (lout,5100) neval, h(i), p
  END IF
END DO

!     START OF MAIN CYCLE.

!     FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).

Main_loop: DO
  loop = loop + 1
  imax = 1
  imin = 1
  hmax = h(1)
  hmin = h(1)
  DO i = 2, np1
    IF (h(i) > hmax) THEN
      imax = i
      hmax = h(i)
    ELSE
      IF (h(i) < hmin) THEN
        imin = i
        hmin = h(i)
      END IF
    END IF
  END DO

!     FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)

  pbar = zero
  DO i = 1, np1
    IF (i /= imax) THEN
      pbar = pbar + g(i,:)
    END IF
  END DO
  pbar = pbar / nap

!     REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
!     HSTAR = FUNCTION VALUE AT PSTAR.

  pstar = a * (pbar - g(imax,:)) + pbar
  CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pstar,hstar)

  neval = neval + 1
  IF (iprint > 0) THEN
    IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, hstar, pstar
  END IF

!     IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
!     HSTST = FUNCTION VALUE AT PSTST.
	
  IF (hstar < hmin) THEN
    pstst = c * (pstar - pbar) + pbar
    CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pstst,hstst)
    neval = neval + 1
    IF (iprint > 0) THEN
      IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, hstst, pstst
    END IF

!     IF HSTST < HMIN REPLACE CURRENT MAXIMUM POINT BY PSTST AND
!     HMAX BY HSTST, THEN TEST FOR CONVERGENCE.

    IF (hstst >= hmin) THEN   ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
      g(imax,:) = pstar
      h(imax) = hstar
    ELSE
      g(imax,:) = pstst
      h(imax) = hstst
    END IF
    GO TO 250
  END IF

!     HSTAR IS NOT < HMIN.
!     TEST WHETHER IT IS < FUNCTION VALUE AT SOME POINT OTHER THAN
!     P(IMAX).   IF IT IS REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.

  DO i = 1, np1
    IF (i /= imax) THEN
      IF (hstar < h(i)) THEN  ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
        g(imax,:) = pstar
        h(imax) = hstar
        GO TO 250
      END IF
    END IF
  END DO

!     HSTAR > ALL FUNCTION VALUES EXCEPT POSSIBLY HMAX.
!     IF HSTAR <= HMAX, REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.

  IF (hstar <= hmax) THEN
    g(imax,:) = pstar
    hmax = hstar
    h(imax) = hstar
  END IF
  
!     CONTRACTED STEP TO THE POINT PSTST,
!     HSTST = FUNCTION VALUE AT PSTST.

  pstst = b * g(imax,:) + (one-b) * pbar
  CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pstst,hstst)
  neval = neval + 1
  IF (iprint > 0) THEN
    IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, hstst, pstst
  END IF

!     IF HSTST < HMAX REPLACE P(IMAX) BY PSTST & HMAX BY HSTST.

  IF (hstst <= hmax) THEN
    g(imax,:) = pstst
    h(imax) = hstst
    GO TO 250
  END IF

!     HSTST > HMAX.
!     SHRINK THE SIMPLEX BY REPLACING EACH POINT, OTHER THAN THE CURRENT
!     MINIMUM, BY A POINT MID-WAY BETWEEN ITS CURRENT POSITION AND THE
!     MINIMUM.

  DO i = 1, np1
    IF (i /= imin) THEN
      DO j = 1, nop
        IF (step(j) /= zero) g(i,j) = (g(i,j) + g(imin,j)) * half
        p(j) = g(i,j)
      END DO
      CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,h(i))
      neval = neval + 1
      IF (iprint > 0) THEN
        IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, h(i), p
      END IF
    END IF
  END DO

!     IF LOOP = NLOOP TEST FOR CONVERGENCE, OTHERWISE REPEAT MAIN CYCLE.

  250 IF (loop < nloop) CYCLE Main_loop

!     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
!     CURRENT SIMPLEX.

  hmean = SUM( h(1:np1) ) / np1
  hstd = SUM( (h(1:np1) - hmean) ** 2 )
  hstd = SQRT(hstd / np1)

!     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
!     START OF THE MAIN CYCLE AGAIN.

  IF (hstd > stopcr .AND. neval <= maxfn) THEN
    iflag = 0
    loop = 0
    CYCLE Main_loop
  END IF

!     FIND THE CENTROID OF THE CURRENT SIMPLEX AND THE FUNCTION VALUE THERE.

  DO i = 1, nop
    IF (step(i) /= zero) THEN
      p(i) = SUM( g(1:np1,i) ) / np1
    END IF
  END DO
  CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,func)
  neval = neval + 1
  IF (iprint > 0) THEN
    IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, func, p
  END IF

!     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, maxfn, HAS BEEN
!     OVERRUN; IF SO, EXIT WITH IFAULT = 1.

  IF (neval > maxfn) THEN
    ifault = 1
    IF (iprint < 0) RETURN
    WRITE (lout,5200) maxfn
    WRITE (lout,5300) hstd
    WRITE (lout,5400) p
    WRITE (lout,5500) func
    RETURN
  END IF

!     CONVERGENCE CRITERION SATISFIED.
!     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
!     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.

  IF (iprint >= 0) THEN
    WRITE (lout,5600)
    WRITE (lout,5400) p
    WRITE (lout,5500) func
  END IF

  IF (iflag == 0 .OR. ABS(savemn-hmean) >= stopcr) THEN
    iflag = 1
    savemn = hmean
    loop = 0
  ELSE
    EXIT Main_loop
  END IF

END DO Main_loop

IF (iprint >= 0) THEN
  WRITE (lout,5700) neval
  WRITE (lout,5800) p
  WRITE (lout,5900) func
END IF

!CAP: modified to skip surface fitting when the minimum is outside
!or within 1% from the edges of the grid
if (iquad <= 0 .or. minval(p)<= 0.01_dp .or. maxval(p)>= 0.99_dp) then
        if (iprint >= 0) then 
          write(lout,*)'minim: solution within 1% from grid edges'
          write(lout,*)'       SKIPPING surface fitting.'
        endif
        RETURN
endif
!or when all the points in the simplex have collapsed onto one spot
do i=1,nop
        if (maxval(g(1:np1,i)-g(1,i)) < tiny(1._dp)) then
                if (iprint >= 0) then 
                  write(lout,*)'minim: points in simplex have collapsed.'
                  write(lout,*)'       SKIPPING surface fitting.' 
                endif
                RETURN
        endif
enddo

!------------------------------------------------------------------

!     QUADRATIC SURFACE FITTING

IF (iprint >= 0) WRITE (lout,6000)

!     EXPAND THE FINAL SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
!     ERRORS.

hmin = func
nmore = 0
DO i = 1, np1
  DO
    test = ABS(h(i)-func)
    IF (test < simp) THEN
      DO j = 1, nop
        IF (step(j) /= zero) g(i,j) = (g(i,j)-p(j)) + g(i,j)
        pstst(j) = g(i,j)
      END DO
      CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pstst,h(i))
      nmore = nmore + 1
      neval = neval + 1

      !CAP: check that we don't exceed the max. number of evaluations
      IF (neval > maxfn) THEN
        if (iprint >= 0) then 
          WRITE (lout,*) 'minim: cannot complete surface fitting'
          WRITE (lout,*) '        neval=',neval,'>',maxfn
        endif
        RETURN
      END IF

      IF (h(i) >= hmin) CYCLE
      hmin = h(i)
      IF (iprint >= 0) WRITE (lout,5100) neval, hmin, pstst
    ELSE
      EXIT
    END IF
  END DO
END DO


!     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.

DO i = 1, nap
  i1 = i + 1
  pstar = (g(1,:) + g(i1,:)) * half
  CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pstar,aval(i))
  nmore = nmore + 1
  neval = neval + 1
END DO


!     THE MATRIX OF ESTIMATED SECOND DERIVATIVES IS CALCULATED AND ITS
!     LOWER TRIANGLE STORED IN BMAT.

a0 = h(1)
DO i = 1, nap
  i1 = i - 1
  i2 = i + 1
  DO j = 1, i1
    j1 = j + 1
    pstst = (g(i2,:) + g(j1,:)) * half
    CALL objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pstst,hstst)
    nmore = nmore + 1
    neval = neval + 1
    l = i * (i-1) / 2 + j
    bmat(l) = two * (hstst + a0 - aval(i) - aval(j))
  END DO
END DO

l = 0
DO i = 1, nap
  i1 = i + 1
  l = l + i
  bmat(l) = two * (h(i1) + a0 - two*aval(i))
END DO

!     THE VECTOR OF ESTIMATED FIRST DERIVATIVES IS CALCULATED AND
!     STORED IN AVAL.

DO i = 1, nap
  i1 = i + 1
  aval(i) = two * aval(i) - (h(i1) + 3._dp*a0) * half
END DO

!     THE MATRIX Q OF NELDER & MEAD IS CALCULATED AND STORED IN G.

pmin = g(1,:)
DO i = 1, nap
  i1 = i + 1
  g(i1,:) = g(i1,:) - g(1,:)
END DO

DO i = 1, nap
  i1 = i + 1
  g(i,:) = g(i1,:)
END DO

!     INVERT BMAT

CALL syminv(bmat, nap, bmat, temp, nullty, ifault, rmax)
IF (ifault == 0) THEN
  irank = nap - nullty
ELSE                                 ! BMAT not +ve definite
                                     ! Resume search for the minimum
  IF (iprint >= 0) WRITE (lout,6100)
  ifault = 2
  IF (neval > maxfn) RETURN
  WRITE (lout,6200)
  step = half * step
  GO TO 20
END IF

!     BMAT*A/2 IS CALCULATED AND STORED IN H.

DO i = 1, nap
  h(i) = zero
  DO j = 1, nap
    IF (j <= i) THEN
      l = i * (i-1) / 2 + j
    ELSE
      l = j * (j-1) / 2 + i
    END IF
    h(i) = h(i) + bmat(l) * aval(j)
  END DO
END DO

!     FIND THE POSITION, PMIN, & VALUE, YMIN, OF THE MINIMUM OF THE
!     QUADRATIC.

ymin = DOT_PRODUCT( h(1:nap), aval(1:nap) )
ymin = a0 - ymin
DO i = 1, nop
  pstst(i) = DOT_PRODUCT( h(1:nap), g(1:nap,i) )
END DO
pmin = pmin - pstst
IF (iprint >= 0) THEN
  WRITE (lout,6300) ymin, pmin
  WRITE (lout,6400)
END IF

!     Q*BMAT*Q'/2 IS CALCULATED & ITS LOWER TRIANGLE STORED IN VC

DO i = 1, nop
  DO j = 1, nap
    h(j) = zero
    DO k = 1, nap
      IF (k <= j) THEN
        l = j * (j-1) / 2 + k
      ELSE
        l = k * (k-1) / 2 + j
      END IF
      h(j) = h(j) + bmat(l) * g(k,i) * half
    END DO
  END DO

  DO j = i, nop
    l = j * (j-1) / 2 + i
    vc(l) = DOT_PRODUCT( h(1:nap), g(1:nap,j) )
  END DO
END DO

!     THE DIAGONAL ELEMENTS OF VC ARE COPIED INTO VAR.

j = 0
DO i = 1, nop
  j = j + i
  var(i) = vc(j)
END DO
IF (iprint < 0) RETURN
WRITE (lout,6500) irank
CALL print_tri_matrix(vc, nop, lout)

WRITE (lout,6600)
CALL syminv(vc, nap, bmat, temp, nullty, ifault, rmax)

!     BMAT NOW CONTAINS THE INFORMATION MATRIX

WRITE (lout,6700)
CALL print_tri_matrix(bmat, nop, lout)

ii = 0
ij = 0
DO i = 1, nop
  ii = ii + i
  IF (vc(ii) > zero) THEN
    vc(ii) = one / SQRT(vc(ii))
  ELSE
    vc(ii) = zero
  END IF
  jj = 0
  DO j = 1, i - 1
    jj = jj + j
    ij = ij + 1
    vc(ij) = vc(ij) * vc(ii) * vc(jj)
  END DO
  ij = ij + 1
END DO

WRITE (lout,6800)
ii = 0
DO i = 1, nop
  ii = ii + i
  IF (vc(ii) /= zero) vc(ii) = one
END DO
CALL print_tri_matrix(vc, nop, lout)

!     Exit, on successful termination.

WRITE (lout,6900) nmore
RETURN

5000 FORMAT (' Progress Report every',i4,' function evaluations'/  &
             ' EVAL.   FUNC.VALUE.          PARAMETER VALUES')
5100 FORMAT (/' ', i4, '  ', g12.5, '  ', 5G11.4, 3(/t22, 5G11.4))
5200 FORMAT (' No. of function evaluations > ',i5)
5300 FORMAT (' RMS of function values of last simplex =', g14.6)
5400 FORMAT (' Centroid of last simplex =',4(/' ', 6G13.5))
5500 FORMAT (' Function value at centroid =', g14.6)
5600 FORMAT (/' EVIDENCE OF CONVERGENCE')
5700 FORMAT (/' Minimum found after',i5,' function evaluations')
5800 FORMAT (' Minimum at',4(/' ', 6G13.6))
5900 FORMAT (' Function value at minimum =', g14.6)
6000 FORMAT (/' Fitting quadratic surface about supposed minimum'/)
6100 FORMAT (/' MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.'/ &
             ' MINIMUM PROBABLY NOT FOUND'/)
6200 FORMAT (/t11, 'Search restarting'/)
6300 FORMAT (' Minimum of quadratic surface =',g14.6,' at',4(/' ', 6G13.5))
6400 FORMAT (' IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED ',      &
             'FROM THE MINIMIZATION,'/' THE MINIMUM MAY BE FALSE &/OR THE ',  &
             'INFORMATION MATRIX MAY BE INACCURATE'/)
6500 FORMAT (' Rank of information matrix =',i3/ &
             ' Inverse of information matrix:-')
6600 FORMAT (/' If the function minimized was -LOG(LIKELIHOOD),'/         &
             ' this is the covariance matrix of the parameters.'/         &
             ' If the function was a sum of squares of residuals,'/       &
             ' this matrix must be multiplied by twice the estimated ',   &
             'residual variance'/' to obtain the covariance matrix.'/)
6700 FORMAT (' INFORMATION MATRIX:-'/)
6800 FORMAT (/' CORRELATION MATRIX:-')
6900 FORMAT (/' A further',i4,' function evaluations have been used'/)

END SUBROUTINE minim




SUBROUTINE syminv(a, n, c, w, nullty, ifault, rmax)

!     ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968.

!     ARGUMENTS:-
!     A()    = INPUT, THE SYMMETRIC MATRIX TO BE INVERTED, STORED IN
!                LOWER TRIANGULAR FORM
!     N      = INPUT, ORDER OF THE MATRIX
!     C()    = OUTPUT, THE INVERSE OF A (A GENERALIZED INVERSE IF C IS
!                SINGULAR), ALSO STORED IN LOWER TRIANGULAR FORM.
!                C AND A MAY OCCUPY THE SAME LOCATIONS.
!     W()    = WORKSPACE, DIMENSION AT LEAST N.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, APPROXIMATE BOUND ON THE ACCURACY OF THE DIAGONAL
!                ELEMENTS OF C.  E.G. IF RMAX = 1.E-04 THEN THE DIAGONAL
!                ELEMENTS OF C WILL BE ACCURATE TO ABOUT 4 DEC. DIGITS.

!     LATEST REVISION - 1 April 1985

!***************************************************************************

REAL (dp), INTENT(IN OUT) :: a(:), c(:), w(:)
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(OUT)      :: nullty, ifault
REAL (dp), INTENT(OUT)    :: rmax

REAL (dp), PARAMETER :: zero = 0._dp, one = 1._dp
INTEGER              :: i, icol, irow, j, jcol, k, l, mdiag, ndiag, nn, nrow
REAL (dp)            :: x

nrow = n
ifault = 1
IF (nrow > 0) THEN
  ifault = 0

!     CHOLESKY FACTORIZATION OF A, RESULT IN C
  CALL chola(a, nrow, c, nullty, ifault, rmax, w)
  IF (ifault == 0) THEN

!     INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERSE
!     OF C, ROW BY ROW STARTING WITH THE LAST ROW.
!     IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW.

    nn = nrow * (nrow+1) / 2
    irow = nrow
    ndiag = nn
    10 IF (c(ndiag) /= zero) THEN
      l = ndiag
      DO i = irow, nrow
        w(i) = c(l)
        l = l + i
      END DO
      icol = nrow
      jcol = nn
      mdiag = nn

      30 l = jcol
      x = zero
      IF (icol == irow) x = one / w(irow)
      k = nrow
      40 IF (k /= irow) THEN
        x = x - w(k) * c(l)
        k = k - 1
        l = l - 1
        IF (l > mdiag) l = l - k + 1
        GO TO 40
      END IF

      c(l) = x / w(irow)
      IF (icol == irow) GO TO 60
      mdiag = mdiag - icol
      icol = icol - 1
      jcol = jcol - 1
      GO TO 30
    END IF ! (c(ndiag) /= zero)

    l = ndiag
    DO j = irow, nrow
      c(l) = zero
      l = l + j
    END DO

    60 ndiag = ndiag - irow
    irow = irow - 1
    IF (irow /= 0) GO TO 10
  END IF
END IF
RETURN

END SUBROUTINE syminv



SUBROUTINE chola(a, n, u, nullty, ifault, rmax, r)

!     ALGORITHM AS6, APPLIED STATISTICS, VOL.17, 1968, WITH
!     MODIFICATIONS BY A.J.MILLER

!     ARGUMENTS:-
!     A()    = INPUT, A +VE DEFINITE MATRIX STORED IN LOWER-TRIANGULAR
!                FORM.
!     N      = INPUT, THE ORDER OF A
!     U()    = OUTPUT, A LOWER TRIANGULAR MATRIX SUCH THAT U*U' = A.
!                A & U MAY OCCUPY THE SAME LOCATIONS.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, AN ESTIMATE OF THE RELATIVE ACCURACY OF THE
!                DIAGONAL ELEMENTS OF U.
!     R()    = OUTPUT, ARRAY CONTAINING BOUNDS ON THE RELATIVE ACCURACY
!                OF EACH DIAGONAL ELEMENT OF U.

!     LATEST REVISION - 1 April 1985

!***************************************************************************

REAL (dp), INTENT(IN)   :: a(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: u(:)
INTEGER, INTENT(OUT)    :: nullty, ifault
REAL (dp), INTENT(OUT)  :: rmax, r(:)

!     ETA SHOULD BE SET EQUAL TO THE SMALLEST +VE VALUE SUCH THAT
!     1._dp + ETA IS CALCULATED AS BEING GREATER THAN 1._dp IN THE ACCURACY
!     BEING USED.

REAL (dp), PARAMETER :: eta = EPSILON(1.0_dp), zero = 0._dp
INTEGER              :: i, icol, irow, j, k, l, m
REAL (dp)            :: rsq, w

ifault = 1
IF (n > 0) THEN
  ifault = 2
  nullty = 0
  rmax = eta
  r(1) = eta
  j = 1
  k = 0

!     FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.

  DO  icol = 1, n
    l = 0

!     IROW = ROW NUMBER WITHIN COLUMN ICOL

    DO  irow = 1, icol
      k = k + 1
      w = a(k)
      IF (irow == icol) rsq = (w*eta) ** 2
      m = j
      DO  i = 1, irow
        l = l + 1
        IF (i == irow) EXIT
        w = w - u(l) * u(m)
        IF (irow == icol) rsq = rsq + (u(l)**2*r(i)) ** 2
        m = m + 1
      END DO

      IF (irow == icol) EXIT
      IF (u(l) /= zero) THEN
        u(k) = w / u(l)
      ELSE
        u(k) = zero
        IF (ABS(w) > ABS(rmax*a(k))) GO TO 60
      END IF
    END DO

!     END OF ROW, ESTIMATE RELATIVE ACCURACY OF DIAGONAL ELEMENT.

    rsq = SQRT(rsq)
    IF (ABS(w) > 5.*rsq) THEN
      IF (w < zero) RETURN
      u(k) = SQRT(w)
      r(i) = rsq / w
      IF (r(i) > rmax) rmax = r(i)
    ELSE
      u(k) = zero
      nullty = nullty + 1
    END IF

    j = j + icol
  END DO
  ifault = zero
END IF
60 RETURN

END SUBROUTINE chola



SUBROUTINE print_tri_matrix(a, n, lout)

INTEGER, INTENT(IN)    :: n, lout
REAL (dp), INTENT(IN)  :: a(:)

!     Local variables
INTEGER  :: i, ii, i1, i2, l

l = 1
DO l = 1, n, 6
  ii = l * (l-1) / 2
  DO i = l, n
    i1 = ii + l
    ii = ii + i
    i2 = MIN(ii,i1+5)
    WRITE (lout,'(1X,6G13.5)') a(i1:i2)
  END DO
END DO
RETURN

END SUBROUTINE print_tri_matrix

end module nelder_mead
