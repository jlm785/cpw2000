C*
      SUBROUTINE SCFFTB (N,C,WSAVE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(1), WSAVE(1)
C
      IF (N .EQ. 1) RETURN
C
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL SCFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
C
      RETURN
      END
