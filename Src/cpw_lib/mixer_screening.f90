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

!>  Calculates a model inverse dielectric function
!>  Should be checked when time allows...
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         22 September 2015, 8 October 2025.
!>  \copyright    GNU Public License v2

  subroutine mixer_screening(xmix, linit, ek, ztot, vcell)

! Written 22 September 2015, extracted from mixer_bfgs. JLM
! Modified, documentation, January 2020. JLM
! Modified, indentation, ztot=0. JLM
! Modified, unphysical values of xmix. 8 October 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  logical, intent(in)                ::  linit                           !<  initialization flag
  real(REAL64), intent(in)           ::  ek                              !<  Kinetic energy (Hartree). G**2 / 2, must be positive
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  vcell                           !<  cell volume

! output

  real(REAL64), intent(out)          ::  xmix                            !<  mixing coefficient

! local variables

  real(REAL64), save        ::  qtf2, alpha, fac, qf                     !  k_TF**2, Slater alpha, xc-screening, k_F
  real(REAL64)              ::  x                                        !  k / k_F
  real(REAL64)              ::  gx                                       !  Lindhard response function
  real(REAL64)              ::  axo2, xl

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter     ::  SMALL = 1.0E-10_REAL64


! ztot=0 in empty lattice

  if(ztot < SMALL) then
    xmix = 0.3
  else
    if(linit) then

      qtf2 = 4.0*(3.0*ztot/(PI*vcell))**(1.0/3.0)
      alpha = 0.77
      fac = -6.0*alpha/(PI*PI*qtf2)
      qf = PI*qtf2/4.0

    else

      x = sqrt(2*ek)/qf
      axo2 = abs(2*UM-x)
      gx = UM / 2
      if (axo2 /= ZERO) then
        xl = log(axo2/(2*UM+x))
        gx = (UM-(UM-x*x/(4*x)) * xl) / 2
      endif
      xmix = UM / (UM + (fac + qtf2/(2*ek)) * gx)

      xmix = abs(xmix)
      if(xmix > UM) xmix = UM

    endif

  endif

  return

end subroutine mixer_screening
