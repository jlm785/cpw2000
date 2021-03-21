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

       SUBROUTINE RFFT_C16(CHD,NN1,N1,N2,N3,MODE,WRK,MXDWRK)
C
C      COMPUTES A REAL TO COMPLEX 3D FAST FOURIER TRANSFORM
C      USING THE CFFT_MLT_C16 INTERFACE SUBROUTINE
C      WHICH IS DIFFERENT FOR FFTW, SCILIB, TROULLIER, NCAR
C      FFT PACKAGES
c
C      INPUT:
C      CHD(2,I,J,K)  3D ARRAY TO BE TRANSFORMED. IT IS COMPLEX,
C                    SO CARE SHOULD BE TAKEN THAT IT IS THE RIGHT FORMAT.
C      NN1           LEADING DIMENSION OF CHD, USUALLY N1+1 (DIRTY TRICK)
C      N1,N2,N3      GRID SIZE < NMFFT (RECOMPILE IF YOU NEED IT TO BE BIGGER)
C      MODE          DETERMINES TRANSFORM MODE: MODE=1 FORWARD, RENORMALIZES
C                    REAL TO COMPLEX (DIRECT-> RECIPROCAL).
C                    MODE=-1 BACKWARD, COMPLEX TO REAL (RECIPROCAL->DIRECT)
C
C      OUTPUT:
C      CHD(2,I,J,K)  TRANSFORMED 3D ARRAY
C
C      WORK:
C      WRK(I)        SIZE MAY DEPEND ON FFT PACKAGE
C                    2*MAX(N1,N2,N3)*MIDDLE(N1,N2,N3) SHOULD BE ENOUGH
C
C      WRITTEN JULY 18, 2002. JLM
C      MODIFIED MARCH 31 2004. JLM
C      MODIFIED 23 FEBRUARY 2006. JLM
C      COPYRIGHT INESC-MN/JOSE LUIS MARTINS
C
C      VERSION 4.44
C
       IMPLICIT NONE
C
C
C      MAXIMUM ALLOWED SIZE IS NMFFT. 2**NMFAX SHOULD BE BIGGER THAN NMFFT.
C
       INTEGER NMFFT,NMFAX
       PARAMETER(NMFFT = 8192, NMFAX = 21)
C
       INTEGER N1,N2,N3,NN1,MODE,MXDWRK
       REAL*8 CHD(2,NN1,N2,N3),WRK(MXDWRK)
C
       INTEGER PLAN2
       INTEGER*8 PLAN3

       INTEGER IFAX(NMFAX)
       INTEGER I,J,K
       INTEGER ISTRIDE,IWORD,JSKIP
       REAL*8 FAC
       REAL*8 TRIGS(2*NMFFT)
       INTEGER N2J2,N3K2,N1I2
       REAL*8 SR,SI,AR,AI
C
       REAL*8 ZERO,UM,DOIS
       PARAMETER (UM = 1.0D0, ZERO=0.0D0, DOIS=2.0D0)
C
       SAVE PLAN2,PLAN3
C
C      INITIAL TESTS
C
       IF(N1 .GT. NMFFT .OR. N2 .GT. NMFFT .OR. N3 .GT. NMFFT) THEN
         WRITE(6,301) N1,N2,N3,NMFFT
         STOP
       ENDIF
       IF(MODE .NE. 1 .AND. MODE .NE. -1) THEN
         WRITE(6,300) MODE
         STOP
       ENDIF
       IF(N1 .LT. 2 .OR. N2 .LT. 2 .OR. N3 .LT. 3) THEN
         WRITE(6,302) N1,N2,N3
         STOP
       ENDIF
       IF(2*MAX(N1*N2,N1*N3,N2*N3) .GT. MXDWRK) THEN
         WRITE(6,303) MXDWRK,2*N1*N2,2*N1*N3,2*N2*N3
         STOP
       ENDIF
C
       IF (MODE .EQ. 1) THEN
C
C        3RD DIM-N3:   PACKS TWO REAL ARRAYS AS A COMPLEX ARRAY
C
         ISTRIDE = NN1*N2
         IWORD = 1
         JSKIP = (N2+1)/2
         FAC = UM/DOIS
C
         CALL CFFT_PREPARE(CHD,WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)

!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,jskip,
!$omp&        plan2,plan3,fac)
         DO J=1,N2/2
           DO K=1,N3
           DO I=1,N1
             CHD(2,I,J,K) = CHD(1,I,J+JSKIP,K)
           ENDDO
           ENDDO
C
           CALL CFFT_MLT_C16(CHD(1,1,J,1),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)
C
C          UNPACKS THE RESULT AND SCALES WHILE IT IS IN CACHE
C
           K=1
           FAC = UM / DBLE(N1*N2*N3)
           DO I=1,N1
             CHD(1,I,J,K) = CHD(1,I,J,K) * FAC
             CHD(1,I,J+JSKIP,K) = CHD(2,I,J,K) * FAC
             CHD(2,I,J,K) = ZERO
             CHD(2,I,J+JSKIP,K) = ZERO
           ENDDO
C
           FAC = UM / DBLE(2*N1*N2*N3)
           IF(N3 .GT. 2) THEN
             DO K=2,(N3+1)/2
               N3K2 = N3 - K + 2
               DO I=1,N1
                 SR = (CHD(1,I,J,K) + CHD(1,I,J,N3K2)) * FAC
                 SI = (CHD(2,I,J,K) + CHD(2,I,J,N3K2)) * FAC
                 AR = (CHD(1,I,J,K) - CHD(1,I,J,N3K2)) * FAC
                 AI = (CHD(2,I,J,K) - CHD(2,I,J,N3K2)) * FAC
                 CHD(1,I,J,K) = SR
                 CHD(2,I,J,K) = AI
C                CHD(1,I,J,N3K2) = SR
C                CHD(2,I,J,N3K2) =-AI
                 CHD(1,I,J+JSKIP,K) = SI
                 CHD(2,I,J+JSKIP,K) =-AR
C                CHD(1,I,J+JSKIP,N3K2) = SI
C                CHD(2,I,J+JSKIP,N3K2) = AR
               ENDDO
             ENDDO
           ENDIF
C
C          NYQUIST
C
           IF(MOD(N3,2) .EQ. 0) THEN
             K=N3/2+1
             FAC = UM / DBLE(N1*N2*N3)
             DO I=1,N1
               CHD(1,I,J,K) = CHD(1,I,J,K) * FAC
               CHD(1,I,J+JSKIP,K) = CHD(2,I,J,K) * FAC
               CHD(2,I,J,K) = ZERO
               CHD(2,I,J+JSKIP,K) = ZERO
             ENDDO
           ENDIF
         ENDDO
!$omp end parallel do
C
C        CORRECTS FOR ODD N2
C
         IF(MOD(N2,2) .EQ. 1) THEN
           DO K=1,N3
           DO I=1,N1
             CHD(2,I,JSKIP,K) = ZERO
           ENDDO
           ENDDO
C
           CALL CFFT_MLT_C16(CHD(1,1,JSKIP,1),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)

           FAC = UM / DBLE(N1*N2*N3)
           DO K=1,N3
           DO I=1,N1
             CHD(1,I,JSKIP,K) = CHD(1,I,JSKIP,K) * FAC
             CHD(2,I,JSKIP,K) = CHD(2,I,JSKIP,K) * FAC
           ENDDO
           ENDDO

         ENDIF

         CALL CFFT_FINISH(PLAN2,PLAN3)
C
C        2ND DIM N2 : DO COMPLEX FFT FOR ALL - NOTE HALF OF N3
C
         ISTRIDE = NN1
         IWORD = 1
         CALL CFFT_PREPARE(CHD,WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N2,N1,MODE,PLAN2,PLAN3)
!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,plan2,plan3)
         DO K=1,N3/2+1
           CALL CFFT_MLT_C16(CHD(1,1,1,K),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N2,N1,MODE,PLAN2,PLAN3)
         ENDDO
!$omp end parallel do

         CALL CFFT_FINISH(PLAN2,PLAN3)
C
C        1ST DIM N1 : DO COMPLEX FFT FOR N1 - NOTE HALF OF N3
C
         ISTRIDE = 1
         IWORD = NN1
         CALL CFFT_PREPARE(CHD,WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)

         K=1
         CALL CFFT_MLT_C16(CHD(1,1,1,K),WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)

!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,plan2,plan3)
         DO K=2,(N3+1)/2
           CALL CFFT_MLT_C16(CHD(1,1,1,K),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)
C
C          UNPACK 001 DIRECTION (WHILE DATA IS IN CACHE)
C
           N3K2 = N3 - K + 2
           CHD(1,1,1,N3K2) =  CHD(1,1,1,K)
           CHD(2,1,1,N3K2) = -CHD(2,1,1,K)
C
C          UNPACK 011 PLANE AND 101 PLANE (MINUS 001 DIRECTIN)
C
           DO J=2,N2
             N2J2 = N2 - J + 2
             CHD(1,1,N2J2,N3K2) =  CHD(1,1,J,K)
             CHD(2,1,N2J2,N3K2) = -CHD(2,1,J,K)
           ENDDO
           DO I=2,N1
             N1I2 = N1 - I + 2
             CHD(1,N1I2,1,N3K2) =  CHD(1,I,1,K)
             CHD(2,N1I2,1,N3K2) = -CHD(2,I,1,K)
           ENDDO
C
C          UNPACK REST
C
           DO J=2,N2
             N2J2 = N2 - J + 2
             DO I=2,N1
               N1I2 = N1 - I + 2
               CHD(1,N1I2,N2J2,N3K2) =  CHD(1,I,J,K)
               CHD(2,N1I2,N2J2,N3K2) = -CHD(2,I,J,K)
             ENDDO
           ENDDO
C
         ENDDO
!$omp end parallel do

         IF(MOD(N3,2) .EQ. 0) THEN
           K = N3/2+1
           CALL CFFT_MLT_C16(CHD(1,1,1,K),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)
         ENDIF
C
         CALL CFFT_FINISH(PLAN2,PLAN3)

       ELSE
C
C        BACKWARD MODE
C
C        1ST DIM - N1 : COMPLEX, ONLY HALF USED
C
         ISTRIDE = 1
         IWORD = NN1
         CALL CFFT_PREPARE(CHD,WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)
!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,plan2,plan3)
         DO K=1,N3/2+1
           CALL CFFT_MLT_C16(CHD(1,1,1,K),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)
         ENDDO
!$omp end parallel do
C

         CALL CFFT_FINISH(PLAN2,PLAN3)
C
C        2ND DIM - N2 : COMPLEX, ONLY HALF USED
C
         ISTRIDE = NN1
         IWORD = 1
         CALL CFFT_PREPARE(CHD,WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N2,N1,MODE,PLAN2,PLAN3)
!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,plan2,plan3)
         DO K=1,N3/2+1
           CALL CFFT_MLT_C16(CHD(1,1,1,K),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N2,N1,MODE,PLAN2,PLAN3)
         ENDDO
!$omp end parallel do
C

         CALL CFFT_FINISH(PLAN2,PLAN3)
C
C        3RD DIM - N3 :  PACKS THE COMPLEX CONJUGATE
C                        AND DOES TWO FFTS AT A TIME
C
         ISTRIDE = NN1*N2
         IWORD = 1
         JSKIP = (N2+1)/2
           CALL CFFT_PREPARE(CHD,WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)
!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,jskip,
!$omp&        plan2,plan3)
         DO J=1,N2/2
           DO K=N3/2+2,N3
           N3K2 = N3 - K + 2
           DO I=1,N1
             CHD(1,I,J,K) = CHD(1,I,J,N3K2) + CHD(2,I,J+JSKIP,N3K2)
             CHD(2,I,J,K) =-CHD(2,I,J,N3K2) + CHD(1,I,J+JSKIP,N3K2)
           ENDDO
           ENDDO
           DO K=1,N3/2+1
           DO I=1,N1
               CHD(1,I,J,K) = CHD(1,I,J,K) - CHD(2,I,J+JSKIP,K)
               CHD(2,I,J,K) = CHD(2,I,J,K) + CHD(1,I,J+JSKIP,K)
           ENDDO
           ENDDO
C
           CALL CFFT_MLT_C16(CHD(1,1,J,1),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)
C
           DO K=1,N3
           DO I=1,N1
             CHD(1,I,J+JSKIP,K) = CHD(2,I,J,K)
             CHD(2,I,J,K) = ZERO
             CHD(2,I,J+JSKIP,K) = ZERO
           ENDDO
           ENDDO

         ENDDO
!$omp end parallel do
C
C        CORRECTS FOR ODD N2
C
         IF(MOD(N2,2) .EQ. 1) THEN
C
           J = JSKIP
           DO K=N3/2+2,N3
           N3K2 = N3 - K + 2
           DO I=1,N1
             CHD(1,I,J,K) = CHD(1,I,J,N3K2)
             CHD(2,I,J,K) =-CHD(2,I,J,N3K2)
           ENDDO
           ENDDO
           DO K=1,N3/2+1
           DO I=1,N1
             CHD(1,I,J,K) = CHD(1,I,J,K)
             CHD(2,I,J,K) = CHD(2,I,J,K)
           ENDDO
           ENDDO

           CALL CFFT_MLT_C16(CHD(1,1,JSKIP,1),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)

           DO K=1,N3
           DO I=1,N1
             CHD(2,I,JSKIP,K) = ZERO
           ENDDO
           ENDDO
         ENDIF
C
         CALL CFFT_FINISH(PLAN2,PLAN3)

       ENDIF
C
       RETURN
 300   FORMAT('  *** STOPPED IN RFFT_C16     DIMENSION FOR FFT= ',
     1         3I5,' EXCEEDS ',I5)
 301   FORMAT('  *** STOPPED IN RFFT_C16     MODE = ',I5,
     2        ' SHOULD BE 1 OR -1')
 302   FORMAT('  *** STOPPED IN RFFT_C16     DIMENSION FOR FFT= ',
     1         3I5,' LESS THEN 2 (OR 3 FOR N3)')
 303   FORMAT('  *** STOPPED IN RFFT_C16     DIMENSION FOR WRKFFT= ',
     1         I5,' LESS THAN MAXIMUM OF',3I5)
       END
