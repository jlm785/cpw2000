       SUBROUTINE CFTFAX (NFFT, IFAX, TRIG)
C
C      WRITTEN JULY 24, 1992 BY NORMAN J. TROULLIER JR.
C      COPYRIGHT NORMAN J. TROULLIER JR.
C      DISTRIBUTION ONLY BY WRITTEN PERMISSION
C
       IMPLICIT NONE
C
       INTEGER IFAX(*),NFFT
C
       REAL*8 TRIG(*)
C
       INTEGER I,II(8), II3, II4, II5, II6
C
       REAL*8 ANG, ANG1, RR1, RR2
C
       LOGICAL LL1, LL2
C
       ANG = 8.D0 * ATAN (1.D0) / DBLE (NFFT)
       DO I=1,NFFT
         ANG1 = ANG * DBLE (I - 1)
         TRIG(I*2-1) = COS (ANG1)
         TRIG(I*2) = SIN (ANG1)
       ENDDO
C
       II(1) = 10
       II(2) = 9
       II(3) = 8
       II(4) = 6
       II(5) = 5
       II(6) = 4
       II(7) = 3
       II(8) = 2
       LL1 = .FALSE.
       LL2 = .FALSE.
       DO I=1,8
         IF (NFFT .EQ. II(I)) LL1 = .TRUE.
       ENDDO
       IF (LL1) THEN
         IFAX(2) = NFFT
         IFAX(1) = 1
       ELSE
         II3 = NFFT
         IFAX(1) = 0
    2    RR1 = 2.D0
    3    RR2 = DBLE (II3) ** (1.D0 / RR1)
         IF (RR2 .GT. DBLE (II(1))) THEN
           RR1 = RR1 + 1.D0
           GOTO 3
         ENDIF
         II4 = INT (RR2)
    5    DO I=1,8
           IF (II4 .EQ. II(I)) GO TO 6
         ENDDO
         II4 = II4 - 1
         IF (II4 .EQ. 1) STOP ' FACTORIZATION FAILED IN CFTFAX'
         GOTO 5
    6    DO II5=I,8
           IF ((II3 / II(II5)) * II(II5) .EQ. II3) THEN
             IFAX(1) = IFAX(1) + 1
             IFAX(IFAX(1)+1) = II(II5)
             II3 = II3 / II(II5)
             DO II6=1,8
               IF (II3 .EQ. II(II6)) LL2 = .TRUE.
             ENDDO
             IF (LL2) THEN
               IFAX(1) = IFAX(1) + 1
               IFAX(IFAX(1)+1) = II3
               GOTO 7
             ENDIF
             GOTO 2
           ENDIF
         ENDDO
         STOP ' FACTORIZATION FAILED IN CFTFAX'
       ENDIF
 7     CONTINUE
       RETURN
       END
