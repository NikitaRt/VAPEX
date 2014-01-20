!---------------------------------------------------------------------!
!  Elliptic solver (BiCGSTAB) for pressure correction                 !
!---------------------------------------------------------------------!
!  $Id: MOD_Elleq.f90 8 2014-01-14 13:38:05Z Sergey $
!---------------------------------------------------------------------!
MODULE ELLEQ
  USE GLOBALS_2D
  !------------------------------------- Multigrid variables
  REAL(KIND=8), ALLOCATABLE:: &
       dpRHS(:,:)
  INTEGER:: NX_ELL, NY_ELL
  !------------------------------------- BiCGSTAB2 variables
  INTEGER,PARAMETER:: L_bcg2=2
  LOGICAL:: okprint_bcg2 
  LOGICAL:: nonzero_bcg2
  CHARACTER(len=3),PARAMETER:: typestop='abs'
  EXTERNAL MATVEC_BCG2
  INTEGER:: info_bcg2
  INTEGER:: mxmv_bcg2
  INTEGER:: ldw_bcg2,n_bcg2
  REAL(KIND=8), ALLOCATABLE:: work_bcg2(:)
  REAL(KIND=8):: tol_bcg2
  !------------------------------------- Boundary conditions
  REAL(KIND=8), ALLOCATABLE:: &
       r1x(:),q1x(:),rnx(:),qnx(:), &
       r1y(:),q1y(:),rmy(:),qmy(:)
  !---------------------------------------------------------
CONTAINS
  SUBROUTINE InitEllipticSolver
    IMPLICIT NONE
    !-------------------------- BiCGSTAB2 solver
    n_bcg2 = n9*m9
    ldw_bcg2 = n_bcg2*(2*l_bcg2+5)
    tol_bcg2 = 1.D-10
    mxmv_bcg2 = 5000
    ALLOCATE(work_bcg2(ldw_bcg2))
    NX_ELL=n9; NY_ELL=m9
    ALLOCATE(r1x(m9),q1x(m9),rnx(m9),qnx(m9), &
         r1y(n9),q1y(n9),rmy(n9),qmy(n9))
    ALLOCATE(dpRHS(n9,m9))
  END SUBROUTINE InitEllipticSolver

  SUBROUTINE bicgstab2 (okprint,l, n, x, rhs, matvec, nonzero, tol, &
       typestop,mxmv, work, ldw, info)

    ! Simple BiCGstab(\ell) iterative method, \ell <= 2
    ! By M.A.Botchev, Jan.'98
    ! report bugs to botchev@cwi.nl or botchev@excite.com

    ! Copyright (c) 1998 by M.A.Botchev
    ! Permission to copy all or part of this work is granted,
    ! provided that the copies are not made or distributed
    ! for resale, and that the copyright notice and this
    ! notice are retained.

    ! This is the "vanilla" version of BiCGstab(\ell) as described
    ! in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements
    ! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
    ! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
    ! properties of BiCGstab methods in finite precision arithmetic",
    ! Numerical Algorithms, 10, 1995, pp.203-223
    ! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
    ! hybrid BiCG methods", Computing, 56, 1996, pp.141-163

    ! {{ This code based on:
    ! subroutine bistbl v1.0 1995

    ! Copyright (c) 1995 by D.R. Fokkema.
    ! Permission to copy all or part of this work is granted,
    ! provided that the copies are not made or distributed
    ! for resale, and that the copyright notice and this
    ! notice are retained.  }}

    ! okprint == (input) LOGICAL. If okprint=.true. residual norm
    ! will be printed to *
    ! l       == (input) INTEGER BiCGstab's dimension <= 2
    ! Set l=2 for highly nonsymmetric problems
    ! n       == (input) INTEGER size of the system to solve
    ! x       == (input/output) REAL(KIND=8) array dimension n
    ! initial guess on input, solution on output
    ! rhs     == (input) REAL(KIND=8) array dimension n
    ! right-hand side (rhs) vector
    ! matvec  == (input) EXTERNAL name of matrix vector subroutine
    ! to deliver y:=A*x by CALL matvec(n,x,y)
    ! nonzero == (input) LOGICAL tells
    ! BiCGstab(\ell) if initial guess in x is zero or not.
    ! If nonzero is .FALSE., one MATVEC call is saved.
    ! tol     == (input/output) REAL(KIND=8) tolerance for all possible
    ! stopping criteria (see the 'typestop' parameter)
    ! On output, if info=0 or 1, tol is actually achieved
    ! residual reduction or whatever (see the 'typestop' parameter)
    ! typestop== (input) CHARACTER*3 stopping criterion (||.|| denotes
    ! the 2-norm):
    ! typestop='rel' -- relative stopping crit.: ||res|| < tol*||res0||
    ! typestop='abs' -- absolute stopping crit.: ||res||<tol
    ! typestop='max' -- maximum  stopping crit.: max(abs(res))<tol
    ! NOTE(for typestop='rel' and 'abs'): To save computational work, the value of
    ! residual norm used to check the convergence inside the main iterative
    ! loop is computed from
    ! projections, i.e. it can be smaller than the true residual norm
    ! (it may happen when e.g. the 'matrix-free' approach is used).
    ! Thus, it is possible that the true residual does NOT satisfy
    ! the stopping criterion ('rel' or 'abs').
    ! The true residual norm (or residual reduction) is reported on
    ! output in parameter TOL -- this can be changed to save 1 MATVEC
    ! (see comments at the end of the subroutine)
    ! mxmv   ==  (input/output) INTEGER.  On input: maximum number of matrix
    ! vector multiplications allowed to be done.  On output:
    ! actual number of matrix vector multiplications done
    ! work   ==  (workspace) REAL(KIND=8) array dimension (n,2*l+5))
    ! ldw    ==  (input) INTEGER size of work, i.e. ldw >= n*(2*l+5)
    ! info   ==  (output) INTEGER.  info = 0 in case of normal computations
    ! and
    ! info = -m (<0) - means paramater number m has an illegal value
    ! info = 1 - means no convergence achieved (stopping criterion
    ! is not fulfilled)
    ! info = 2 - means breakdown of the algorithm (try to enlarge
    ! parameter l=\ell to get rid of this)
    ! ----------------------------------------------------------
    IMPLICIT NONE
    EXTERNAL matvec
    INTEGER :: l, n, mxmv, ldw, info
    REAL(KIND=8)::  x(n), rhs(n), tol
    LOGICAL ::   okprint,nonzero
    CHARACTER(3) ::   typestop
    REAL(KIND=8)::  work(n,5+2*l)
    ! -----------------------------------------
    INTEGER ::   lmax
    PARAMETER(lmax=2)
    REAL(KIND=8):: rwork(lmax+1,3+2*(lmax+1))

    LOGICAL :: GoOn, rcmp, xpdt
    INTEGER :: ii, i1, jj, kk, nmv
    REAL(KIND=8):: alpha,beta,hatgamma,kappa0, kappal,maxval1, &
         mxnrmr,mxnrmx,omega,rho0,rho1,rnrm0,rnrm,rnrmMax, &
         sigma,sum1,varrho
    INTEGER :: z, zz, y0, yl, y
    INTEGER :: rr, r, u, xp, bp

    REAL(KIND=8)::    zero, one, delta
    PARAMETER(zero=0d0,one=1d0,delta=1d-2)

    REAL(KIND=8),PARAMETER:: Tiny = 1.d-25
    REAL(KIND=8)::  rho_1_0

    info = 0

    IF (l > lmax .OR. l < 1) info = -2
    IF (tol <= zero) info = -8
    IF (mxmv < 0) info = -10

    rr = 1
    r = rr+1
    u = r+(l+1)
    xp = u+(l+1)
    bp = xp+1
    IF (bp*n > ldw) info = -12

    z = 1
    zz = z+(l+1)
    y0 = zz+(l+1)
    yl = y0+1
    y = yl+1

    IF (info /= 0) RETURN

    ! Initialize first residual
    IF (nonzero) THEN
       CALL matvec (n, x, work(1,r) )
       DO ii=1,n
          work(ii,r) = rhs(ii) - work(ii,r)
       ENDDO
       nmv = 1
    ELSE
       DO ii=1,n
          work(ii,r) = rhs(ii)
       ENDDO
       nmv = 0
    ENDIF

    ! Initialize iteration loop
    sum1 = zero
    DO ii=1,n
       work(ii,rr) = work(ii,r)
       work(ii,bp) = work(ii,r)
       work(ii,xp) = x(ii)
       x(ii) = zero
       sum1=sum1+ work(ii,r)**2
    ENDDO

    rnrm0 = SQRT( sum1 )
    rnrm = rnrm0
    IF (typestop == 'max') THEN
       maxval1=zero
       DO ii=1,n
          maxval1=MAX( maxval1, ABS( work(ii,r) ) )
       ENDDO
       rnrmMax = maxval1
    ENDIF

    mxnrmx = rnrm0
    mxnrmr = rnrm0
    rcmp = .FALSE.
    xpdt = .FALSE.

    alpha = zero
    omega = one
    sigma = one
    rho0 =  one

    ! Iterate
    IF(typestop == 'rel')THEN
       GoOn = rnrm.GE.tol*rnrm0 .AND. nmv.LT.mxmv
    ELSE IF(typestop == 'abs')THEN
       GoOn = rnrm.GE.tol       .AND. nmv.LT.mxmv
    ELSE IF(typestop == 'max')THEN
       GoOn = rnrmMax.GE.tol    .AND. nmv.LT.mxmv
    ELSE
       info = -9
       RETURN
    END IF

    DO WHILE (GoOn)
       ! =====================
       ! --- The BiCG part ---
       ! =====================
       rho0 = -omega*rho0
       DO kk=1,l
          sum1 = zero
          DO ii=1,n
             sum1=sum1+ work(ii,rr)*work(ii,r+kk-1)
          ENDDO
          rho1 = sum1

          !          if (rho0 == zero) then
          !             info = 2
          !             return
          !          endif

          IF (ABS(rho0).LT.Tiny) THEN
             !     SY: check, maybe we got solution
             DO ii=1,n
                x(ii) = work(ii,xp) + x(ii)
             ENDDO

             CALL matvec (n, x, work(1,r) )
             sum1 = zero
             DO ii=1,n
                work(ii,r)=work(ii,r)-rhs(ii)
                sum1=sum1+work(ii,r)**2
             ENDDO
             rnrm = SQRT( ABS(sum1) )

             IF(typestop.EQ.'rel'.AND.rnrm.LT.tol*rnrm0) GOTO 999
             IF(typestop.EQ.'abs'.AND.rnrm.LT.tol) GOTO 999
             IF(typestop.EQ.'max')THEN
                maxval1=zero
                DO ii=1,n
                   maxval1=MAX(maxval1, ABS( work(ii,r) ) )
                ENDDO
                rnrmMax = maxval1
                IF(rnrmMax.LT.tol) GOTO 999
             ENDIF

             info = 2
             RETURN
          ENDIF

          rho_1_0 = MIN(rho1/rho0,10.d1)
          rho_1_0 = MAX(rho_1_0,-10.d1)

          beta = alpha*rho_1_0 !(rho1/rho0)


          !          beta = alpha*(rho1/rho0)
          rho0 = rho1
          DO jj=0,kk-1
             DO ii=1,n
                work(ii,u+jj) = work(ii,r+jj) - beta*work(ii,u+jj)
             ENDDO
          ENDDO

          CALL matvec(n, work(1,u+kk-1), work(1,u+kk))
          nmv = nmv+1

          sum1 = zero
          DO ii=1,n
             sum1=sum1+ work(ii,rr)*work(ii,u+kk)
          ENDDO
          sigma = sum1

          !          if (sigma == zero) then
          !             info = 2
          !             return
          !          endif
          IF (ABS(sigma).LT.Tiny) THEN

             !     SY: check, maybe we got solution
             DO ii=1,n
                x(ii) = work(ii,xp) + x(ii)
             ENDDO

             CALL matvec (n, x, work(1,r) )
             sum1 = zero
             DO ii=1,n
                work(ii,r)=work(ii,r)-rhs(ii)
                sum1=sum1+work(ii,r)**2
             ENDDO
             rnrm = SQRT( ABS(sum1) )

             IF(typestop.EQ.'rel'.AND.rnrm.LT.tol*rnrm0) GOTO 999
             IF(typestop.EQ.'abs'.AND.rnrm.LT.tol) GOTO 999
             IF(typestop.EQ.'max')THEN
                maxval1=zero
                DO ii=1,n
                   maxval1=MAX(maxval1, ABS( work(ii,r) ) )
                ENDDO
                rnrmMax = maxval1
                IF(rnrmMax.LT.tol) GOTO 999
             ENDIF
             info = 2
             RETURN
          ENDIF

          !      	    alpha = rho1/sigma

          alpha = rho1/sigma
          DO ii=1,n
             x(ii) = alpha*work(ii,u) + x(ii)
          ENDDO

          DO jj=0,kk-1
             DO ii=1,n
                work(ii,r+jj) = -alpha*work(ii,u+jj+1) + work(ii,r+jj)
             ENDDO
          ENDDO

          CALL matvec (n, work(1,r+kk-1), work(1,r+kk))
          nmv = nmv+1

          sum1 = zero
          DO ii=1,n
             sum1=sum1+ work(ii,r)**2
          ENDDO
          rnrm = SQRT( sum1 )
          IF (typestop == 'max') THEN
             maxval1=zero
             DO ii=1,n
                maxval1=MAX( maxval1, ABS( work(ii,r) ) )
             ENDDO
             rnrmMax = maxval1
          ENDIF
          mxnrmx = MAX (mxnrmx, rnrm)
          mxnrmr = MAX (mxnrmr, rnrm)
       ENDDO

       ! ==================================
       ! --- The convex polynomial part ---
       ! ==================================

       ! Z = R'R
       DO i1=1,l+1
          DO jj=i1-1,l
             sum1 = zero
             DO ii=1,n
                sum1=sum1+ work(ii,r+jj)*work(ii,r+i1-1)
             ENDDO
             rwork(jj+1,z+i1-1) = sum1
             rwork(z+i1-1,jj+1) = rwork(jj+1,z+i1-1)
          ENDDO
       ENDDO

       DO i1=zz,zz+l
          DO ii=1,l+1
             rwork(ii,i1)   = rwork(ii,i1+(z-zz))
          ENDDO
       ENDDO
       ! tilde r0 and tilde rl (small vectors)

       rwork(1,y0) = -one
       rwork(2,y0) = rwork(2,z) / rwork(2,zz+1)
       rwork(l+1,y0) = zero

       rwork(1,yl) = zero
       rwork(2,yl) = rwork(2,z+l) / rwork(2,zz+1)
       rwork(l+1,yl) = -one

       ! Convex combination
       DO ii=1,l+1
          rwork(ii,y) = zero
       ENDDO
       DO jj=1,l+1
          DO ii=1,l+1
             rwork(ii,y) = rwork(ii,y) + rwork(jj,yl)* &
                  rwork(ii,z+jj-1)
          ENDDO
       ENDDO
       sum1 = zero
       DO ii=1,l+1
          sum1=sum1+ rwork(ii,yl)*rwork(ii,y)
       ENDDO
       kappal = SQRT( sum1 )

       DO ii=1,l+1
          rwork(ii,y) = zero
       ENDDO
       DO jj=1,l+1
          DO ii=1,l+1
             rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)* &
                  rwork(ii,z+jj-1)
          ENDDO
       ENDDO
       sum1 = zero
       DO ii=1,l+1
          sum1=sum1+ rwork(ii,y0)*rwork(ii,y)
       ENDDO
       kappa0 = SQRT( sum1 )

       sum1 = zero
       DO ii=1,l+1
          sum1=sum1+ rwork(ii,yl)*rwork(ii,y)
       ENDDO
       varrho = sum1
       varrho = varrho / (kappa0*kappal)

       hatgamma = SIGN(1d0,varrho)*MAX(ABS(varrho),7d-1) * &
            (kappa0/kappal)

       DO ii=1,l+1
          rwork(ii,y0) = -hatgamma*rwork(ii,yl) +      rwork(ii,y0)
       ENDDO

       ! Update
       omega = rwork(l+1,y0)
       DO jj=1,l
          DO ii=1,n
             work(ii,u) = work(ii,u) - rwork(1+jj,y0)*work(ii,u+jj)
             x(ii)      = x(ii)      + rwork(1+jj,y0)*work(ii,r+jj-1)
             work(ii,r) = work(ii,r) - rwork(1+jj,y0)*work(ii,r+jj)
          ENDDO
       ENDDO

       DO ii=1,l+1
          rwork(ii,y) = zero
       ENDDO
       DO jj=1,l+1
          DO ii=1,l+1
             rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)* &
                  rwork(ii,z+jj-1)
          ENDDO
       ENDDO

       sum1 = zero
       DO ii=1,l+1
          sum1=sum1+ rwork(ii,y0)*rwork(ii,y)
       ENDDO
       rnrm = SQRT( sum1 )

       ! ================================
       ! --- The reliable update part ---
       ! ================================
       mxnrmx = MAX (mxnrmx, rnrm)
       mxnrmr = MAX (mxnrmr, rnrm)
       xpdt = (rnrm.LT.delta*rnrm0.AND.rnrm0.LT.mxnrmx)
       rcmp = ((rnrm.LT.delta*mxnrmr.AND.rnrm0.LT.mxnrmr) .OR.xpdt)
       IF (rcmp) THEN
          CALL matvec (n, x, work(1,r) )
          nmv = nmv + 1
          DO ii=1,n
             work(ii,r) = work(ii,bp) - work(ii,r)
          ENDDO
          mxnrmr = rnrm
          IF (xpdt) THEN
             DO ii=1,n
                work(ii,xp) = x(ii) + work(ii,xp)
                x(ii) = zero
                work(ii,bp) = work(ii,r)
             ENDDO
             mxnrmx = rnrm
          ENDIF
       ENDIF

       IF(typestop == 'rel')THEN
          GoOn = rnrm.GE.tol*rnrm0 .AND. nmv.LT.mxmv
          IF(okprint) PRINT *, nmv,' ',rnrm/rnrm0
       ELSE IF(typestop == 'abs')THEN
          GoOn = rnrm.GE.tol       .AND. nmv.LT.mxmv
          IF(okprint) PRINT *, nmv,' ',rnrm
       ELSE IF(typestop == 'max')THEN
          maxval1=zero
          DO ii=1,n
             maxval1=MAX(maxval1, ABS( work(ii,r) ) )
          ENDDO
          rnrmMax = maxval1
          GoOn = rnrmMax.GE.tol    .AND. nmv.LT.mxmv
          IF(okprint) PRINT *, nmv,' ',rnrmMax
       END IF

    ENDDO

    ! =========================
    ! --- End of iterations ---
    ! =========================

    DO ii=1,n
       x(ii) = work(ii,xp) + x(ii)
    ENDDO

999 CONTINUE

    ! --------------------- One matvec can be saved by commenting out this:
    ! (this is to compute the true residual)
    CALL matvec (n, x, work(1,r) )
    DO ii=1,n
       work(ii,r) = rhs(ii) - work(ii,r)
    ENDDO

    sum1 = zero
    DO ii=1,n
       sum1=sum1+ work(ii,r)**2
    ENDDO
    rnrm = SQRT( sum1 )
    ! --------------------- One matvec can be saved by commenting out this^

    IF(typestop == 'rel')THEN
       IF (rnrm > tol*rnrm0) info = 1
       tol = rnrm/rnrm0
    ELSE IF(typestop == 'abs')THEN
       IF (rnrm > tol) info = 1
       tol = rnrm
    ELSE IF(typestop == 'max')THEN
       IF (rnrmMax > tol) info = 1
       tol = rnrmMax
    END IF

    mxmv = nmv

    RETURN
  END SUBROUTINE bicgstab2

END MODULE ELLEQ
