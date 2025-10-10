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

!>  Gets the filename for pseudopotential files from path,
!>  chemical symbol and suffix
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         10 October 2025.
!>  \copyright    GNU Public License v2

subroutine read_pseudo_get_path(nameatom, fnam, pseudo_path, pseudo_suffix)

! Written 10 October 2025. JLM


  implicit none

! input

  character(len=2), intent(in)       ::  nameatom                        !<  chemical symbol

  character(len=200), intent(in)     ::  pseudo_path                     !<  path to pseudopotentials
  character(len=50), intent(in)      ::  pseudo_suffix                   !<  suffix for the pseudopotentials

! output

  character(len=255), intent(out)    ::  fnam                            !<  name of file



  if(nameatom(1:1) /= ' ' .and. nameatom(2:2) /= ' ') then
    fnam = nameatom//adjustl(trim(pseudo_suffix))
  else if(nameatom(1:1) == ' ') then
    fnam = nameatom(2:2)//adjustl(trim(pseudo_suffix))
  else
    fnam = nameatom(1:1)//adjustl(trim(pseudo_suffix))
  endif
  fnam = adjustl(trim(pseudo_path))//fnam

  return

end subroutine read_pseudo_get_path
