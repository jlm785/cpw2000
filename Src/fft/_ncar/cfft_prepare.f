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

       SUBROUTINE CFFT_PREPARE(A,C,TRIGS,IFAX,ISKIP,ISKIP2,NFFT,
     1 NUMBER,ISIGN,PLAN2,PLAN3)
C
C      INTERFACE SUBROUTINE FOR NCAR
C
C      UNUSED: EVERYTHING!
C
C
C      WRITTEN FEBRUARY 22,2006.JLM
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
c
       CALL SCFFTI(NFFT,TRIGS)

       RETURN
       END
