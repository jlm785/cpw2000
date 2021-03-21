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

       SUBROUTINE CFFT_MLT_C16(A,C,TRIGS,IFAX,ISKIP,ISKIP2,NFFT,  &
            NUMBER,ISIGN,PLAN2,PLAN3)
!
!      INTERFACE SUBROUTINE FOR INTEL MKL LIBRARY DFTI MODULES
!      ONLY A AND ISIGN ARE USED, EVERYTHING ELSE HAS BEEN SET PREVIOUSLY...
!
!      INPUT:
!      NFFT        TRANSFORM LENGTH
!      NUMBER      NUMBER OF FFTS TO BE PERFORMED IN PARALLEL
!      A           COMPLEX ARRAY CONTAINING INPUT
!      C           ADDITIONAL WORK ARRAY OF LENGTH 2 * N * NUMBER
!      TRIGS       COMPLEX ARRAY OF DFT EXPONENTIAL FACTORS
!                  OF LENGTH N (PREVIOUSLY GENERATED)
!      IFAX(1)     NUMBER OF FACTORS IN FFT NUMBER NFFT ACCORDING
!                  TO NFFT=N1*N2....*NK WHERE NK=2,3,4,5
!      IFAX(2)     IFAX(NK+1) - THE FFT FACTORS PREVIOUSLY
!                  CALCULATED (IFAX AND TRIGS CAN BE
!                  CALCULATED BY: CALL CFTFAX(NFFT,IFAX,TRIGS)
!      ISKIP       STRIDE OF FFT WITHIN ARRAYS A AND B
!      ISKIP2      INCREMENT IN WORDS BETWEEN CONSEQUETIVE
!                  TRANSFORMS
!      ISIGN        DETERMINES TRANSFORM MODE MODE=1 FORWARD, RENORMALIZES
!                  (DIRECT-> RECIPRO!AL). MODE=-1 BACKWARD (RECIPROCAL->DIRECT)
!
!      OUTPUT:
!      A           COMPLEX ARRAY CONTAINING OUPUT
!
!      WORK:
!      C           ADDITIONAL WORK ARRAY OF LENGTH 2*NFFT
!
!      WRITTEN MARCH 16, 2010.JLM
!      COPYRIGHT INESC-MN/JOSE LUIS MARTINS
!
!      VERSION 4.50
!
       USE MKL_DFTI
       USE MKL_HANDLE
       IMPLICIT NONE

       integer, parameter          :: REAL64 = selected_real_kind(12)

       INTEGER IFAX(*),ISIGN,ISKIP,ISKIP2,NFFT,NUMBER
       REAL(REAL64)   :: C(2,NFFT),TRIGS(*)
       COMPLEX(REAL64) :: A(*)

       INTEGER PLAN2
       INTEGER*8 PLAN3

       INTEGER STATUS

       IF(ISIGN .EQ. 1) THEN
         Status = DftiComputeForward(hand, A)
       ELSE
         Status = DftiComputeBackward(hand, A)
       ENDIF

       RETURN
       END
