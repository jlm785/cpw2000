       SUBROUTINE CFFTMLT(A,B,C,TRIGS,IFAX,ISKIP,ISKIP2,NFFT,
     1  NUMBER,ISIGN)
C
C       MIXED RADIX FFT ROUTINE BASED ON TEMPERTON ,JOUR.COMP.
C       PHYS. VOL 52 NUMBER 1 OCT 1983 PAGE 1.
C       CODE BY NORMAN J. TROULLIER JR.
C       NOT FOR SALE, ONLY TO BE USED FOR SCIENTIFIC WORK  
C       REPLACEMENT CODE FOR CRAY ROUTINES, NOTE IT DOES
C       NOT DO THE DATA ROTATION THAT THE CRAY ROUTINE DOES,
C       HENCE ONLY HALF THE WORK SPACE(ARRAY C) IS NEEDED.
C       WRITTEN JULY 24, 1992 BY NORMAN J. TROULLIER JR.
C       COPYRIGHT NORMAN J. TROULLIER JR.
C       DISTRIBUTION ONLY BY WRITTEN PERMISSION
C
C       FFTS ARE PERFORMED IN PARALLEL
C
C       VARIABLES:
C           NFFT    -     TRANSFORM LENGTH
C           NUMBER  - NUMBER OF FFTS TO BE PERFORMED IN PARALLEL
C           A       - ARRAY CONTAINING REAL PART OF INPUT/OUTPUT
C           B       - ARRAY CONTAINING IMAGINARY PART OF INPUT/OUTPUT
C           C       - ADDITIONAL WORK ARRAY OF LENGTH 2 * N * NUMBER
C           TRIGS   - COMPLEX ARRAY OF DFT EXPONENTIAL FACTORS
C                     OF LENGTH N (PREVIOUSLY GENERATED)
C           IFAX(1) - NUMBER OF FACTORS IN FFT NUMBER NFFT ACCORDING
C                     TO NFFT=N1*N2....*NK WHERE NK=2,3,4,5
C           IFAX(2)...IFAX(NK+1) - THE FFT FACTORS PREVIOUSLY
C                     CALCULATED (IFAX AND TRIGS CAN BE 
C                     CALCULATED BY: CALL CFTFAX(NFFT,IFAX,TRIGS)
C           ISKIP  -  STRIDE OF FFT WITHIN ARRAYS A AND B
C           ISKIP2 - INCREMENT IN WORDS BETWEEN CONSEQUETIVE
C                          TRANSFORMS
C           ISIGN   - FFT SIGN
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION A(*),B(*),C(*),IFAX(*),TRIGS(*)
C
       ISKIP3 = NUMBER + NUMBER
       LA = 1
       II = 1
       NFAC = IFAX(1)
       IF(NFAC .NE. 1) THEN
         DO I=1,NFAC
           IFAC = IFAX(I+1)
           IF (II .GT. 0) THEN
             CALL VPASS1(A,B,C,C(1+NUMBER),TRIGS,TRIGS(2),IFAC,LA,
     1        NFFT,ISKIP,ISKIP3,ISKIP2,NUMBER,ISIGN)
           ELSE
             CALL VPASS2(C,C(1+NUMBER),A,B,TRIGS,TRIGS(2),IFAC,LA,
     1        NFFT,ISKIP3,ISKIP,ISKIP2,NUMBER,ISIGN)
           ENDIF
           II = -II
           LA = LA*IFAC
         ENDDO
C
         IF (II .LT. 0) THEN 
           N2 = NFFT*ISKIP
           JJ = 1-ISKIP3
           DO J=1,N2,ISKIP
             JJ=JJ+ISKIP3
             IJ = 0
             JI = JJ+NUMBER
             DO I=0,NUMBER-1
               A(J+IJ) = C(JJ+I)
               B(J+IJ) = C(JI+I)
               IJ = IJ + ISKIP2
             ENDDO
           ENDDO
         ENDIF
       ELSE
         CALL VPASSI(A,B,NFFT,ISKIP,ISKIP2,NUMBER,ISIGN)
       ENDIF
       RETURN
       END