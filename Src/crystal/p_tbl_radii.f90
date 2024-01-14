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

!>  Gives the Clementi atomic radius of an element
!>  JChemPhys 47, 1300 (1967) until radon, above that
!>  educated guess
!>
!>  \author       Jose Luis Martins
!>  \version      5.10
!>  \date         21 April 2021, 29 December 2023.
!>  \copyright    GNU Public License v2


subroutine p_tbl_radii(name, atradius)

! Written  21 April 2021. JLM
! Indentation, ZZ element. 29 December 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len=2), intent(in)       ::  name                            !< chemical symbol

! output

  real(REAL64), intent(out)          ::  atradius                        !< atomic mass (in a.u.)



  if (name == 'H ' .or. name == ' H') then
    atradius = 0.53
  elseif (name == 'D ' .or. name == ' D') then
    atradius = 0.53
  elseif (name == 'T ' .or. name == ' T') then
    atradius = 0.53
  elseif (name == 'He') then
    atradius = 0.31
  elseif (name == 'Li') then
    atradius = 1.67
  elseif (name == 'Be') then
    atradius = 1.12
  elseif (name == 'B ' .or. name == ' B') then
    atradius = 0.87
  elseif (name == 'C ' .or. name == ' C') then
    atradius = 0.67
  elseif (name == 'N ' .or. name == ' N') then
    atradius = 0.56
  elseif (name == 'O ' .or. name == ' O') then
    atradius = 0.48
  elseif (name == 'F ' .or. name == ' F') then
    atradius = 0.42
  elseif (name == 'Ne') then
    atradius = 0.38
  elseif (name == 'Na') then
    atradius = 1.90
  elseif (name == 'Mg') then
    atradius = 1.45
  elseif (name == 'Al') then
    atradius = 1.18
  elseif (name == 'Si') then
    atradius = 1.11
  elseif (name == 'P ' .or. name == ' P') then
    atradius = 0.98
  elseif (name == 'S ' .or. name == ' S') then
    atradius = 0.88
  elseif (name == 'Cl') then
    atradius = 0.79
  elseif (name == 'Ar') then
    atradius = 0.71
  elseif (name == 'K ' .or. name == ' K') then
    atradius = 2.43
  elseif (name == 'Ca') then
    atradius = 1.94
  elseif (name == 'Sc') then
    atradius = 1.84
  elseif (name == 'Ti') then
    atradius = 1.76
  elseif (name == 'V ' .or. name == ' V') then
    atradius = 1.71
  elseif (name == 'Cr') then
    atradius = 1.66
  elseif (name == 'Mn') then
    atradius = 1.61
  elseif (name == 'Fe') then
    atradius = 1.56
  elseif (name == 'Co') then
    atradius = 1.52
  elseif (name == 'Ni') then
    atradius = 1.49
  elseif (name == 'Cu') then
    atradius = 1.45
  elseif (name == 'Zn') then
    atradius = 1.42
  elseif (name == 'Ga') then
    atradius = 1.36
  elseif (name == 'Ge') then
    atradius = 1.25
  elseif (name == 'As') then
    atradius = 1.14
  elseif (name == 'Se') then
    atradius = 1.03
  elseif (name == 'Br') then
    atradius = 0.94
  elseif (name == 'Kr') then
    atradius = 0.88
  elseif (name == 'Rb') then
    atradius = 2.65
  elseif (name == 'Sr') then
    atradius = 2.19
  elseif (name == 'Y ' .or. name == ' Y') then
    atradius = 2.12
  elseif (name == 'Zr') then
    atradius = 2.06
  elseif (name == 'Nb') then
    atradius = 1.98
  elseif (name == 'Mo') then
    atradius = 1.90
  elseif (name == 'Tc') then
    atradius = 1.83
  elseif (name == 'Ru') then
    atradius = 1.78
  elseif (name == 'Rh') then
    atradius = 1.73
  elseif (name == 'Pd') then
    atradius = 1.69
  elseif (name == 'Ag') then
    atradius = 1.65
  elseif (name == 'Cd') then
    atradius = 1.61
  elseif (name == 'In') then
    atradius = 1.56
  elseif (name == 'Sn') then
    atradius = 1.45
  elseif (name == 'Sb') then
    atradius = 1.33
  elseif (name == 'Te') then
    atradius = 1.23
  elseif (name == 'I ' .or. name == ' I') then
    atradius = 1.15
  elseif (name == 'Xe') then
    atradius = 1.08
  elseif (name == 'Cs') then
    atradius = 2.98
  elseif (name == 'Ba') then
    atradius = 2.53
  elseif (name == 'La') then
    atradius = 2.26
  elseif (name == 'Ce') then
    atradius = 2.10
  elseif (name == 'Pr') then
    atradius = 2.47
  elseif (name == 'Nd') then
    atradius = 2.06
  elseif (name == 'Pm') then
    atradius = 2.05
  elseif (name == 'Sm') then
    atradius = 2.38
  elseif (name == 'Eu') then
    atradius = 2.31
  elseif (name == 'Gd') then
    atradius = 2.33
  elseif (name == 'Tb') then
    atradius = 2.25
  elseif (name == 'Dy') then
    atradius = 2.28
  elseif (name == 'Ho') then
    atradius = 2.26
  elseif (name == 'Er') then
    atradius = 2.26
  elseif (name == 'Tm') then
    atradius = 2.22
  elseif (name == 'Yb') then
    atradius = 2.22
  elseif (name == 'Lu') then
    atradius = 2.17
  elseif (name == 'Hf') then
    atradius = 2.08
  elseif (name == 'Ta') then
    atradius = 2.00
  elseif (name == 'W ' .or. name == ' W') then
    atradius = 1.93
  elseif (name == 'Re') then
    atradius = 1.88
  elseif (name == 'Os') then
    atradius = 1.85
  elseif (name == 'Ir') then
    atradius = 1.80
  elseif (name == 'Pt') then
    atradius = 1.77
  elseif (name == 'Au') then
    atradius = 1.74
  elseif (name == 'Hg') then
    atradius = 1.71
  elseif (name == 'Tl') then
    atradius = 1.56
  elseif (name == 'Pb') then
    atradius = 1.54
  elseif (name == 'Bi') then
    atradius = 1.43
  elseif (name == 'Po') then
    atradius = 1.35
  elseif (name == 'At') then
    atradius = 1.27
  elseif (name == 'Rn') then
    atradius = 1.20
  elseif (name == 'Fr') then
    atradius = 3.6
  elseif (name == 'Ra') then
    atradius = 3.0
  elseif (name == 'Ac') then
    atradius = 2.8
  elseif (name == 'Th') then
    atradius = 2.6
  elseif (name == 'Pa') then
    atradius = 2.4
  elseif (name == ' U' .or. name == 'U ') then
    atradius = 2.1
  elseif (name == 'Np') then
    atradius = 2.1
  elseif (name == 'Pu') then
    atradius = 2.3
  elseif (name == 'Am') then
    atradius = 2.3
  elseif (name == 'Cm') then
    atradius = 2.3
  elseif (name == 'Bk') then
    atradius = 2.3
  elseif (name == 'Cf') then
    atradius = 2.3
  elseif (name == 'Es') then
    atradius = 2.3
  elseif (name == 'Fm') then
    atradius = 2.3
  elseif (name == 'Md') then
    atradius = 2.3
  elseif (name == 'No') then
    atradius = 2.3
  elseif (name == 'Lr') then
    atradius = 2.3
  elseif (name == 'ZZ') then
    atradius = 1.0
  else
    atradius = 2.0
  endif

  return

end subroutine p_tbl_radii
