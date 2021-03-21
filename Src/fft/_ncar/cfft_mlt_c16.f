!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

       SUBROUTINE CFFT_MLT_C16(A,C,TRIGS,IFAX,ISKIP,ISKIP2,NFFT,
     1 NUMBER,ISIGN,PLAN2,PLAN3)
C
C      INTERFACE SUBROUTINE FOR CFFTMLT
C
C      INPUT:
C      NFFT        TRANSFORM LENGTH
C      NUMBER      NUMBER OF FFTS TO BE PERFORMED IN PARALLEL
C      A           COMPLEX ARRAY CONTAINING INPUT
C      C           ADDITIONAL WORK ARRAY OF LENGTH 2 * N * NUMBER
C      TRIGS       COMPLEX ARRAY OF DFT EXPONENTIAL FACTORS
C                  OF LENGTH N (PREVIOUSLY GENERATED)
C      IFAX(1)     NUMBER OF FACTORS IN FFT NUMBER NFFT ACCORDING
C                  TO NFFT=N1*N2....*NK WHERE NK=2,3,4,5
C      IFAX(2)     IFAX(NK+1) - THE FFT FACTORS PREVIOUSLY
C                  CALCULATED (IFAX AND TRIGS CAN BE
C                  CALCULATED BY: CALL CFTFAX(NFFT,IFAX,TRIGS)
C      ISKIP       STRIDE OF FFT WITHIN ARRAYS A AND B
C      ISKIP2      INCREMENT IN WORDS BETWEEN CONSEQUETIVE
C                  TRANSFORMS
C      ISIGN        DETERMINES TRANSFORM MODE MODE=1 FORWARD, RENORMALIZES
C                  (DIRECT-> RECIPROCAL). MODE=-1 BACKWARD (RECIPROCAL->DIRECT)
C
C      OUTPUT:
C      A           COMPLEX ARRAY CONTAINING OUPUT
C
C      WORK:
C      C           ADDITIONAL WORK ARRAY OF LENGTH 2*NUMBER*NFFT
C
C      UNUSED PLAN2,PLAN3
C
C      WRITTEN JULY 18, 2002.JLM
C      COPYRIGHT INESC-MN/JOSE LUIS MARTINS
C
C      VERSION 4.40
C
       IMPLICIT NONE
C
       REAL*8 A(2,*),C(*),TRIGS(*)
       INTEGER IFAX(*)
       INTEGER ISKIP,ISKIP2,NFFT,NUMBER,ISIGN

       INTEGER PLAN2
       INTEGER*8 PLAN3
C
       CALL CFFTMLT(A(1,1),A(2,1),C,TRIGS,IFAX,2*ISKIP,2*ISKIP2,NFFT,
     1  NUMBER,ISIGN)


       RETURN
       END
