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

       SUBROUTINE CFFT_PREPARE(A,C,TRIGS,IFAX,ISKIP,ISKIP2,NFFT,   &
                NUMBER,ISIGN,PLAN2,PLAN3)
!
!      INTERFACE SUBROUTINE FOR INTEL MKL LIBRARY DFTI MODULES
!      TRIGS,IFAX,C NOT USED
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
!                  (DIRECT-> RECIPROCAL). MODE=-1 BACKWARD (RECIPROCAL->DIRECT)
!
!      OUTPUT:
!      PLAN3       FFTW3 PLAN
!
!      WORK:
!      C           ADDITIONAL WORK ARRAY OF LENGTH 2*NFFT
!
!      WRITTEN MARCH 16,2010.JLM
!      COPYRIGHT INESC-MN/JOSE LUIS MARTINS
!
!      VERSION 4.52
!
       USE MKL_DFTI
       USE MKL_HANDLE
       IMPLICIT NONE

       INTEGER IFAX(*),ISIGN,ISKIP,ISKIP2,NFFT,NUMBER
       REAL*8 A(2,*),C(2,NFFT),TRIGS(*)

       INTEGER PLAN2
       INTEGER*8 PLAN3
       INTEGER STATUS
       integer   strides_in(2)
       integer   maxthr
       integer   omp_get_max_threads
       external omp_get_max_threads

       strides_in(1) = 0
       strides_in(2) = ISKIP

       maxthr = omp_get_max_threads()

!
!      Create DFTI descriptor
!
       Status = DftiCreateDescriptor(hand, DFTI_DOUBLE, &
                                    DFTI_COMPLEX, 1, NFFT)
!
       Status = DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, NUMBER)
       Status = DftiSetValue(hand, DFTI_INPUT_STRIDES, strides_in)
       Status = DftiSetValue(hand, DFTI_INPUT_DISTANCE, ISKIP2)
       Status = DftiSetValue(hand, DFTI_NUMBER_OF_USER_THREADS, maxthr)
!
!      Commit DFTI descriptor
!
       Status = DftiCommitDescriptor(hand)

       RETURN
       END
