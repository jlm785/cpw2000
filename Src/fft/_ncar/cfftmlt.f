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

        SUBROUTINE CFFTMLT(A,B,C,TRIGS,IFAX,ISKIP,ISKIP2,NFFT,
     1  NUMBER,ISIGN)
C
C       REPLACEMENT CODE FOR CRAY ROUTINES,
C       CALLS NCAR SUBROUTINES WRITTEN BY PAUL N SWARZTRAUBER
C       AND WHICH ARE AVAILABLE FROM NETLIB
C
C       VARIABLES:
C           NFFT    -     TRANSFORM LENGTH
C           NUMBER  - NUMBER OF FFTS TO BE PERFORMED IN PARALLEL
C           A       - ARRAY CONTAINING REAL PART OF INPUT/OUTPUT
C           B       - ARRAY CONTAINING IMAGINARY PART OF INPUT/OUTPUT
C           C       - ADDITIONAL WORK ARRAY
C           TRIGS   - WORK ARRAY OF LENGTH 2*NFFT
C           IFAX    - NOT USED
C           ISKIP   - STRIDE OF FFT WITHIN ARRAYS A AND B
C           ISKIP2  - INCREMENT IN WORDS BETWEEN CONSEQUETIVE
C                          TRANSFORMS
C           ISIGN   - FFT SIGN
C
C      WRITTEN FEBRUARY 10, 2000. JLM
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION A(*),B(*),C(*),IFAX(*),TRIGS(*)
C
       DO I=1,NUMBER
         DO J=1,NFFT
           C(2*J-1) = A((I-1)*ISKIP2+(J-1)*ISKIP+1)
           C(2*J)   = B((I-1)*ISKIP2+(J-1)*ISKIP+1)
         ENDDO
C
         IF(ISIGN .GT. 0) THEN
           CALL SCFFTF(NFFT,C,TRIGS)
         ELSE
           CALL SCFFTB(NFFT,C,TRIGS)
         ENDIF
C
         DO J=1,NFFT
           A((I-1)*ISKIP2+(J-1)*ISKIP+1) = C(2*J-1) 
           B((I-1)*ISKIP2+(J-1)*ISKIP+1) = C(2*J) 
         ENDDO
       ENDDO
C
       RETURN
       END
