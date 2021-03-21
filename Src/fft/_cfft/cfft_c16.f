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

       SUBROUTINE CFFT_C16(CHD,NN1,N1,N2,N3,MODE,WRK,MXDWRK)
C
C      COMPUTES A COMPLEX 3D FAST FOURIER TRANSFORM
C      CALLING THE CFFT_MLT_C16 AND CFTFAX MULTIPLE FFT INTERFACE SUBROUTINE
C
C
C      INPUT:
C      CHD(M,I,J,K)   I=1,N1,J=1,N2,K=1,N3.
C                     M=1 REAL, M=2 IMAGINARY
C                     ARRAY TO BE FOURIER TRANSFORMED
C      NN1            LEADING (SECOND) DIMENSION OF CHD
C                     IN SOME COMPUTERS SHOULD BE ODD (N1+1)
C      N1,N2,N3       GRID SIZE IN THE THREE DIMENSIONS < NMFFT
C      MODE           DETERMINES TRANSFORM MODE
C                     MODE = 1   FORWARD (REAL TO RECIPROCAL) RENORMALIZES
C                     MODE =-1   BACKWARD (RECIPROCAL TO REAL)
C
C      OUTPUT:
C      CHD(M,I,J,K)   I=1,N1,J=1,N2,K=1,N3.
C                     M=1 REAL, M=2 IMAGINARY
C                     ARRAY WITH FOURIER TRANSFORS
C      WORK:
C      WRK(I)         SIZE MAY DEPEND ON FFT PACKAGE
C                     2*MAX(N1,N2,N3)*MIDDLE(N1,N2,N3) SHOULD BE ENOUGH
C
C
C      WRITTEN 26 JULY 2002. JLM
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
C
       INTEGER IFAX(NMFAX)
       INTEGER I,J,K
       INTEGER ISTRIDE,IWORD
       REAL*8 FAC
       REAL*8 TRIGS(2*NMFFT)
C
       REAL*8 UM
       PARAMETER (UM = 1.0D0)

       SAVE PLAN2,PLAN3
C
C
       IF(MODE .NE. 1 .AND. MODE .NE. -1) THEN
         WRITE(6,301) MODE
         STOP
       ENDIF
       IF(N1 .GT. NMFFT .OR. N2 .GT. NMFFT .OR. N3 .GT. NMFFT) THEN
         WRITE(6,300) N1,N2,N3,NMFFT
         STOP
       ENDIF
       IF(N1 .LT. 2 .OR. N2 .LT. 2 .OR. N3 .LT. 2) THEN
         WRITE(6,302) N1,N2,N3
         STOP
       ENDIF
       IF(2*MAX(N1*N2,N1*N3,N2*N3) .GT. MXDWRK) THEN
         WRITE(6,303) MXDWRK,2*N1*N2,2*N1*N3,2*N2*N3
         STOP
       ENDIF
C
       FAC = UM / REAL(N1*N2*N3)

       ISTRIDE = 1
       IWORD = NN1
       CALL CFFT_PREPARE(CHD,WRK,
     1 TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)

       IF(MODE .EQ. 1) THEN
!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,
!$omp&        plan2,plan3,fac)
         DO I=1,N3
           CALL CFFT_MLT_C16(CHD(1,1,1,I),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)
           DO J=1,N2
           DO K=1,N1
             CHD(1,K,J,I) = FAC*CHD(1,K,J,I)
             CHD(2,K,J,I) = FAC*CHD(2,K,J,I)
           ENDDO
           ENDDO
         ENDDO
!$omp end parallel do

       ELSE
!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,plan2,plan3)
         DO I=1,N3
           CALL CFFT_MLT_C16(CHD(1,1,1,I),WRK,
     1     TRIGS,IFAX,ISTRIDE,IWORD,N1,N2,MODE,PLAN2,PLAN3)
         ENDDO
!$omp end parallel do

       ENDIF

       CALL CFFT_FINISH(PLAN2,PLAN3)

C
       ISTRIDE = NN1
       IWORD = 1
       CALL CFFT_PREPARE(CHD,WRK,
     1 TRIGS,IFAX,ISTRIDE,IWORD,N2,N1,MODE,PLAN2,PLAN3)

!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,plan2,plan3)
       DO I=1,N3
         CALL CFFT_MLT_C16(CHD(1,1,1,I),WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N2,N1,MODE,PLAN2,PLAN3)
       ENDDO
!$omp end parallel do

       CALL CFFT_FINISH(PLAN2,PLAN3)

C
       ISTRIDE = NN1*N2
       IWORD = 1
       CALL CFFT_PREPARE(CHD,WRK,
     1 TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)
!$omp parallel do schedule(static,4) default(private)
!$omp& shared(chd,trigs,ifax,istride,iword,n1,n2,n3,mode,plan2,plan3)
       DO I=1,N2
         CALL CFFT_MLT_C16(CHD(1,1,I,1),WRK,
     1   TRIGS,IFAX,ISTRIDE,IWORD,N3,N1,MODE,PLAN2,PLAN3)
       ENDDO
!$omp end parallel do
C
       CALL CFFT_FINISH(PLAN2,PLAN3)
C
       RETURN
 300   FORMAT('  *** STOPPED IN CFFT_C16    DIMENSION FOR FFT= ',
     1         3I5,' EXCEEDS ',I5,' RECOMPILE WITH LARGER DIMENSIONS')
 301   FORMAT('  *** STOPPED IN CFFT_C16     MODE = ',I5,
     2        ' SHOULD BE 1 OR -1')
 302   FORMAT('  *** STOPPED IN CFFT_C16     DIMENSION FOR FFT= ',
     1         3I5,' LESS THAN 2')
 303   FORMAT('  *** STOPPED IN CFFT_C16     DIMENSION FOR WRKFFT= ',
     1         I5,' LESS THAN MAXIMUM OF',3I5)
       END
