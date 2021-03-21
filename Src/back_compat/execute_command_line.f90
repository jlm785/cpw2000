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

!>  compatibility with compilers before fortran 2008,
!>  reverts to call system

subroutine execute_command_line(command)

  implicit none

! input

  character(len=*), intent(in)     ::  command                           !<  command to be executed by system

  call system(command)

end subroutine execute_command_line
