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

!>  Gets the time of day (hh:mm:ss).
!>  Dummy version.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         21 October 2003, 11 September 2015.
!>  \copyright    GNU Public License v2

subroutine zetime(btime)

! Written 21 october 2003
! Modified 11 September 2015. f90. JLM
! Indentation. 28 September 2026. JLM+claude


  implicit none

! output

  character(len=8), intent(out)      ::  btime                           !<  dummy time of day (hh:mm:ss)

  btime = 'NOW     '

  return

end subroutine zetime
