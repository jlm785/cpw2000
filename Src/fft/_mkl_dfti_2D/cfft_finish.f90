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

       SUBROUTINE CFFT_FINISH(PLAN2,PLAN3)
!
!      INTERFACE SUBROUTINE FOR MKL_DFTI SUBROUTINES
!      CLEARS THE 
!
!      UNUSED: PLAN2, PLAN3
!      DATA IS TRANSFERED THROUGH HAND OF MKL_HANDLE MODULE
!
!
!      WRITTEN MARCH 16 2010
!      COPYRIGHT INESC-MN/JOSE LUIS MARTINS
!
!      VERSION 4.50
!
       USE MKL_DFTI
       USE MKL_HANDLE     
       IMPLICIT NONE
!
       INTEGER PLAN2
       INTEGER*8 PLAN3
       INTEGER STATUS

       Status = DftiFreeDescriptor(hand) 
 
       RETURN
       END
