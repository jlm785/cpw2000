       SUBROUTINE LBFGS(N,M,X,G,E,DXMAX,GTOL,CNTROL,H,IH)
c
C Directs a limited memory BFGS minimization of a function which
C is evaluated by the calling program.
c
c The minimization program is from 
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c and this subroutine is a wrapper written 3 October 2006 
c by J.L. Martins to make it the most SIMILAR possible to 
c conjgr subroutine of SIESTA.
c
c  THE DIMENSIONS AND VARIABLES ARE NOT THE SAME!!!!!!!!!!!!!
c
C  N     : INPUT SPACE DIMENSIONALITY
C  M     ; MAXIMUM NUMBER OF VARIABLE METRIC CORRECTIONS
C          USED TO DEFINE THE LIMITED MEMORY MATRIX
C  X     : INPUT POSITION AT WHICH GRADIENT HAS BEEN EVALUATED AND
C          OUTPUT NEW POSITION AT WHICH GRADIENT MUST BE EVALUATED NEXT
C  G     : INPUT GRADIENT (WITH A MINUS SIGN, OR ACTIVATE A LINE BELOW)
C  E     ; ENERGY/ENTHALPY
C  DXMAX : INPUT MAXIMUM ALLOWED DISPLACEMENT
C  GTOL  : INPUT MAXIMUM FINAL VALUE FOR GRADIENT COMPONENT
C  cntrol: MUST BE MADE ZERO BEFORE
C          FIRST CALL. IF IT IS ZERO ON OUTPUT, MINIMIZATION IS
C          CONVERGED. OTHERWISE, CALCULATE GRADIENT AT THE NEW
C          POSITION AND CALL AGAIN THIS ROUTINE. 
C          IF IT IS -1 AN ERROR HAS OCURRED
C
C          DO NOT MODIFY ANY
C          ARGUMENT OTHER THAN G BETWEEN CALLS.
C
C  H,IH  : AUXILIARY ARRAY WHICH MUST NOT BE MODIFIED BETWEEN CALLS
C

       IMPLICIT NONE
C
C      Passed variables
C      
       INTEGER N,M
       REAL*8 X(N),G(N),E,DXMAX(N),GTOL(N),CNTROL
C
C      Passed work variables (should not be changed)
C       
       REAL*8 H( (2*M+7)*N + 11*M*M+8*M )
       INTEGER IH(4*N) 
C
C      Local variables
c       
       CHARACTER*60 TASK,CSAVE
       INTEGER IPRINT,ISAVE(44)
       LOGICAL LSAVE(4)
       REAL*8 DSAVE(29)
       REAL*8 FACTR,GTX
       INTEGER I,ICONV
       
       SAVE CSAVE,DSAVE,ISAVE,LSAVE
       SAVE TASK,IPRINT,FACTR
C
       IF (NINT(CNTROL) .EQ. 0) THEN
         DO I=1,N
           H(N+I) = 0.0D0
           H(2*N+I) = 0.0D0
           IH(I) = 0
         ENDDO
         FACTR = 1.0D7
         TASK = 'START'
         IPRINT = -1        
       ENDIF
C
C      CHANGES THE METRIC SO THAT TRANSFORMED DXMAX IS IDENTITY
C       
       GTX = GTOL(1)*DXMAX(1)
       DO I=1,N
         X(I) = X(I)/DXMAX(I)
         H(I) = -G(I)*DXMAX(I)
         IF(GTOL(I)*DXMAX(I) .LT. GTX) GTX = GTOL(I)*DXMAX(I)
       ENDDO
       
       
       DO I=1,10
         IF(TASK(1:2) .EQ. 'FG') GOTO 503
       
          CALL SETULB(N,M,X,H(N+1),H(2*N+1),IH(1),E,H(1),FACTR,
     1                GTX,H(3*N+1),IH(N+1),TASK,IPRINT,
     2                CSAVE,LSAVE,ISAVE,DSAVE)
     
       ENDDO
       
 503   CONTINUE
       
          CALL SETULB(N,M,X,H(N+1),H(2*N+1),IH(1),E,H(1),FACTR,
     1                GTX,H(3*N+1),IH(N+1),TASK,IPRINT,
     2                CSAVE,LSAVE,ISAVE,DSAVE)
     
       DO I=1,10
         IF(TASK(1:5) .NE. 'NEW_X' ) GOTO 504
       
          CALL SETULB(N,M,X,H(N+1),H(2*N+1),IH(1),E,H(1),FACTR,
     1                GTX,H(3*N+1),IH(N+1),TASK,IPRINT,
     2                CSAVE,LSAVE,ISAVE,DSAVE)
            
       ENDDO
       
 504   CONTINUE
 
       DO I=1,N
         X(I) = X(I)*DXMAX(I)
       ENDDO
       
       IF(TASK(1:2) .EQ. 'FG' .OR. TASK(1:5) .EQ. 'NEW_X') THEN
         CNTROL = 1.0
         ICONV = 1
         DO I=1,N
           IF(ABS(G(I)) .GT. GTOL(I)) ICONV = 0
         ENDDO
         IF(ICONV .EQ. 1) CNTROL = 0.0
       ELSEIF(TASK(1:4) .EQ. 'CONV') THEN
         CNTROL = 0.0
       ELSE
         CNTROL =-1.0
       ENDIF

       RETURN
       END
