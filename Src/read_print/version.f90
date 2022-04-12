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

!>  hardcoded library version
!>
!>  \author       Jose Luis Martins
!>  \version      5.02
!>  \date         13 September 2021.
!>  \copyright    GNU Public License v2

subroutine version(cpwversion, ldevel)

!  Written October 12, 2018. jlm
!  Modified February 2019. JLM
!  Copyright INESC-MN/Jose Luis Martins

   implicit none

!  output

   character(len=4), intent(out)      ::  cpwversion                 !<  hardcoded library version
   logical, intent(out)               ::  ldevel                     !<  development branch, minor version may be incompatible

   cpwversion = '5.05'
   ldevel = .FALSE.

   return

end subroutine version

