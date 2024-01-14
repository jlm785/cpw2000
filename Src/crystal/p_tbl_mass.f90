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

!>  Gives the atomic mass of an element
!>
!>  \author       Jose Luis Martins
!>  \version      5.10
!>  \date         23 October 93, 29 December 2023.
!>  \copyright    GNU Public License v2

subroutine p_tbl_mass(name,atmass)

! Written 23 October 93. JLM
! Modified for f90, 8 June 2014. JLM
! Indentation, ZZ element. 29 December 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len=2), intent(in)       ::  name                            !< chemical symbol

! output

  real(REAL64), intent(out)          ::  atmass                          !< atomic mass (in a.u.)



  if (name == 'H ' .or. name == ' H') then
    atmass = 1.00079_REAL64
  elseif (name == 'D ' .or. name == ' D') then
    atmass = 2.0_REAL64
  elseif (name == 'T ' .or. name == ' T') then
    atmass = 3.0_REAL64
  elseif (name == 'He') then
    atmass = 4.00260_REAL64
  elseif (name == 'Li') then
    atmass = 6.941_REAL64
  elseif (name == 'Be') then
    atmass = 9.01218_REAL64
  elseif (name == 'B ' .or. name == ' B') then
    atmass = 10.81_REAL64
  elseif (name == 'C ' .or. name == ' C') then
    atmass = 12.011_REAL64
  elseif (name == 'N ' .or. name == ' N') then
    atmass = 14.0067_REAL64
  elseif (name == 'O ' .or. name == ' O') then
    atmass = 15.9994_REAL64
  elseif (name == 'F ' .or. name == ' F') then
    atmass = 18.998403_REAL64
  elseif (name == 'Ne') then
    atmass = 20.179_REAL64
  elseif (name == 'Na') then
    atmass = 22.98977_REAL64
  elseif (name == 'Mg') then
    atmass = 24.305_REAL64
  elseif (name == 'Al') then
    atmass = 26.98154_REAL64
  elseif (name == 'Si') then
    atmass = 28.0855_REAL64
  elseif (name == 'P ' .or. name == ' P') then
    atmass = 30.97376_REAL64
  elseif (name == 'S ' .or. name == ' S') then
    atmass = 32.06_REAL64
  elseif (name == 'Cl') then
    atmass = 35.453_REAL64
  elseif (name == 'Ar') then
    atmass = 39.948_REAL64
  elseif (name == 'K ' .or. name == ' K') then
    atmass = 39.0983_REAL64
  elseif (name == 'Ca') then
    atmass = 40.08_REAL64
  elseif (name == 'Sc') then
    atmass = 44.9559_REAL64
  elseif (name == 'Ti') then
    atmass = 47.90_REAL64
  elseif (name == 'V ' .or. name == ' V') then
    atmass = 50.9415_REAL64
  elseif (name == 'Cr') then
    atmass = 51.996_REAL64
  elseif (name == 'Mn') then
    atmass = 54.9380_REAL64
  elseif (name == 'Fe') then
    atmass = 55.847_REAL64
  elseif (name == 'Co') then
    atmass = 58.9332_REAL64
  elseif (name == 'Ni') then
    atmass = 58.70_REAL64
  elseif (name == 'Cu') then
    atmass = 63.546_REAL64
  elseif (name == 'Zn') then
    atmass = 65.38_REAL64
  elseif (name == 'Ga') then
    atmass = 69.72_REAL64
  elseif (name == 'Ge') then
    atmass = 72.59_REAL64
  elseif (name == 'As') then
    atmass = 74.9216_REAL64
  elseif (name == 'Se') then
    atmass = 78.96_REAL64
  elseif (name == 'Br') then
    atmass = 79.904_REAL64
  elseif (name == 'Kr') then
    atmass = 83.80_REAL64
  elseif (name == 'Rb') then
    atmass = 85.4678_REAL64
  elseif (name == 'Sr') then
    atmass = 87.62_REAL64
  elseif (name == 'Y ' .or. name == ' Y') then
    atmass = 88.9059_REAL64
  elseif (name == 'Zr') then
    atmass = 91.22_REAL64
  elseif (name == 'Nb') then
    atmass = 92.9064_REAL64
  elseif (name == 'Mo') then
    atmass = 95.94_REAL64
  elseif (name == 'Tc') then
    atmass = 98.0_REAL64
  elseif (name == 'Ru') then
    atmass = 101.07_REAL64
  elseif (name == 'Rh') then
    atmass = 102.9055_REAL64
  elseif (name == 'Pd') then
    atmass = 106.4_REAL64
  elseif (name == 'Ag') then
    atmass = 107.868_REAL64
  elseif (name == 'Cd') then
    atmass = 112.41_REAL64
  elseif (name == 'In') then
    atmass = 114.82_REAL64
  elseif (name == 'Sn') then
    atmass = 118.69_REAL64
  elseif (name == 'Sb') then
    atmass = 121.75_REAL64
  elseif (name == 'Te') then
    atmass = 127.60_REAL64
  elseif (name == 'I ' .or. name == ' I') then
    atmass = 126.9045_REAL64
  elseif (name == 'Xe') then
    atmass = 131.30_REAL64
  elseif (name == 'Cs') then
    atmass = 132.9054_REAL64
  elseif (name == 'Ba') then
    atmass = 137.33_REAL64
  elseif (name == 'La') then
    atmass = 138.9055_REAL64
  elseif (name == 'Ce') then
    atmass = 140.12_REAL64
  elseif (name == 'Pr') then
    atmass = 140.9077_REAL64
  elseif (name == 'Nd') then
    atmass = 144.24_REAL64
  elseif (name == 'Pm') then
    atmass = 145.0_REAL64
  elseif (name == 'Sm') then
    atmass = 150.4_REAL64
  elseif (name == 'Eu') then
    atmass = 151.96_REAL64
  elseif (name == 'Gd') then
    atmass = 157.25_REAL64
  elseif (name == 'Tb') then
    atmass = 158.9254_REAL64
  elseif (name == 'Dy') then
    atmass = 162.50_REAL64
  elseif (name == 'Ho') then
    atmass = 164.9304_REAL64
  elseif (name == 'Er') then
    atmass = 167.26_REAL64
  elseif (name == 'Tm') then
    atmass = 168.9342_REAL64
  elseif (name == 'Yb') then
    atmass = 173.04_REAL64
  elseif (name == 'Lu') then
    atmass = 174.967_REAL64
  elseif (name == 'Hf') then
    atmass = 178.49_REAL64
  elseif (name == 'Ta') then
    atmass = 180.9479_REAL64
  elseif (name == 'W ' .or. name == ' W') then
    atmass = 183.85_REAL64
  elseif (name == 'Re') then
    atmass = 186.207_REAL64
  elseif (name == 'Os') then
    atmass = 190.2_REAL64
  elseif (name == 'Ir') then
    atmass = 192.22_REAL64
  elseif (name == 'Pt') then
    atmass = 195.09_REAL64
  elseif (name == 'Au') then
    atmass = 196.9665_REAL64
  elseif (name == 'Hg') then
    atmass =200.59_REAL64
  elseif (name == 'Tl') then
    atmass = 204.37_REAL64
  elseif (name == 'Pb') then
    atmass = 207.2_REAL64
  elseif (name == 'Bi') then
    atmass = 208.9804_REAL64
  elseif (name == 'Po') then
    atmass = 209.0_REAL64
  elseif (name == 'At') then
    atmass = 210.0_REAL64
  elseif (name == 'Rn') then
    atmass = 222.0_REAL64
  elseif (name == 'Fr') then
    atmass = 223.0_REAL64
  elseif (name == 'Ra') then
    atmass = 226.0254_REAL64
  elseif (name == 'Ac') then
    atmass = 227.0278_REAL64
  elseif (name == 'Th') then
    atmass = 232.0381_REAL64
  elseif (name == 'Pa') then
    atmass = 231.0359_REAL64
  elseif (name == ' U' .or. name == 'U ') then
    atmass = 238.029_REAL64
  elseif (name == 'Np') then
    atmass = 237.0482_REAL64
  elseif (name == 'Pu') then
    atmass = 244.0_REAL64
  elseif (name == 'Am') then
    atmass = 243.0_REAL64
  elseif (name == 'Cm') then
    atmass = 247.0_REAL64
  elseif (name == 'Bk') then
    atmass = 247.0_REAL64
  elseif (name == 'Cf') then
    atmass = 251.0_REAL64
  elseif (name == 'Es') then
    atmass = 252.0_REAL64
  elseif (name == 'Fm') then
    atmass = 257.0_REAL64
  elseif (name == 'Md') then
    atmass = 258.0_REAL64
  elseif (name == 'No') then
    atmass = 259.0_REAL64
  elseif (name == 'Lr') then
    atmass = 260.0_REAL64
  elseif (name == 'ZZ') then
    atmass = 1000.0_REAL64
  else
    write(6,'("  STOPPED in p_tbl_mass, element ",a2," unknown")') name

    stop

  endif

  return

end subroutine p_tbl_mass
