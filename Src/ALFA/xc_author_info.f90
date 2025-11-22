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

!>  Sets flags according to the exchange and correlation functionals.
!>  Concentrates dispersed conditional statements for ease
!>  of adding new functionals.
!>
!>  \author       José Luís Martins
!>  \version      5.12
!>  \date         22 November 2025.
!>  \copyright    GNU Public License v2

subroutine chr_author_info(author, lxcgrad, lxclap, lxctau,              &
       lxctb09, lxccalc)

! Written August 5, 2002. jlm

  implicit none

! input

  character(len=4), intent(in)       ::  author                          !<  type of xc wanted (ca=pz , pw92 , pbe,...)

!output

  logical, intent(out)               ::  lxcgrad                         !<  gradient of charge density should be calculated
  logical, intent(out)               ::  lxclap                          !<  laplacian of charge density should be calculated
  logical, intent(out)               ::  lxctau                          !<  Kinetic energy density should be calculated
  logical, intent(out)               ::  lxctb09                         !<  Tran-Blaha constant is present
  logical, intent(out)               ::  lxccalc                         !<  xc energy is calculted

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

  lxcgrad = .FALSE.
  lxclap = .FALSE.
  lxctau = .FALSE.
  lxctb09 = .FALSE.
  lxccalc = .TRUE.

  if(chrsameinfo(author, 'PBE') .or. chrsameinfo(author, 'LAK') .or.     &
     chrsameinfo(author, 'TBL') .or. chrsameinfo(author, 'TB09')) then
       lxcgrad = .TRUE.
  endif

  if(chrsameinfo(author, 'TBL') .or. chrsameinfo(author, 'TB09')) then
       lxclap = .TRUE.
       lxctb09 = .TRUE.
       lxccalc = .FALSE.
  endif

  if(chrsameinfo(author, 'TBL') .or. chrsameinfo(author, 'TB09') .or.    &
     chrsameinfo(author, 'LAK')) then
       lxctau = .TRUE.
  endif

  return

end subroutine chr_author_info
