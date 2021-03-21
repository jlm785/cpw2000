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
C      INTERFACE SUBROUTINE FOR VFFT
C
C      INPUT:
C      NFFT        TRANSFORM LENGTH
C
C      OUTPUT:
C      PLAN3       TRIGS,IFAX
C
C      UNUSED: A,C,ISKIP,ISKIP2,NUMBER,ISIGN,PLAN2,PLAN3
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

       CALL CFTFAX (NFFT, IFAX, TRIGS)

       RETURN
       END
