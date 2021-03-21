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
C      INTERFACE SUBROUTINE FOR FFTW3 (WWW.FFTW.ORG)
C      TRIGS,IFAX,C NOT USED
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
C      C           ADDITIONAL WORK ARRAY OF LENGTH 2*NFFT
C
C      WRITTEN DECEMBER 8, 2003.JLM
C      MODIFIED APRIL 10, 2004. JLM
C      MODIFIED FEBRUARY 22, 2006. JLM
C      COPYRIGHT INESC-MN/JOSE LUIS MARTINS
C
C      VERSION 4.43
C
       IMPLICIT NONE
C
       INTEGER IFAX(*),ISIGN,ISKIP,ISKIP2,NFFT,NUMBER
       REAL*8 A(2,*),C(2,NFFT),TRIGS(*)
C
       INTEGER PLAN2
       INTEGER*8 PLAN3

       CALL DFFTW_EXECUTE_DFT(PLAN3,A,A)

       RETURN
       END
