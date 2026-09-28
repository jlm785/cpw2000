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

!>  Interface subroutine for the Intel MKL library DFTI modules.
!>  Frees the descriptor stored in the mkl_handle module.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         16 March 2010.
!>  \copyright    GNU Public License v2

SUBROUTINE CFFT_FINISH(PLAN2, PLAN3)
!
! INTERFACE SUBROUTINE FOR MKL_DFTI SUBROUTINES
! CLEARS THE
!
! UNUSED: PLAN2, PLAN3
! DATA IS TRANSFERED THROUGH HAND OF MKL_HANDLE MODULE
!
!
! WRITTEN MARCH 16 2010
! INDENTATION. 28 SEPTEMBER 2026. JLM+CLAUDE
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

END SUBROUTINE CFFT_FINISH
