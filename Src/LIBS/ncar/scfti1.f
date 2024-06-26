C*
      SUBROUTINE SCFTI1 (N,WA,IFAC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION WA(1)
      DIMENSION IFAC(*), NTRYH(4)
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4) /3, 4, 2, 5/
      DATA TPI   /  6.2831853071 7958647692 5286766559 00577D0/
C
      NL = N
      NF = 0
      J = 0
C
  101 J = J+1
      IF (J.LE.4) NTRY = NTRYH(J)
      IF (J.GT.4) NTRY = NTRY + 2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR.NE.0) GO TO 101
C
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
C
  107 IF (NL .NE. 1) GO TO 104
C
      IFAC(1) = N
      IFAC(2) = NF
C
      ARGH = TPI/REAL(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
C
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.D0
            WA(I) = 0.D0
            LD = LD+L1
            FI = 0.D0
            ARGLD = REAL(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
C
         L1 = L2
  110 CONTINUE
C
      RETURN
      END
