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

!>  Kinetic energy of a plane-wave modified near the cutoff so
!>  that the total energy is a smooth function of the cell shape,
!>  and its derivative with respect to the unmodified kinetic energy
!>  (for the stress, which then includes the Pulay stress).
!>
!>  \author       claude
!>  \version      5.13
!>  \date         15 October 2026.
!>  \copyright    GNU Public License v2


!   The idea is from M. Bernasconi, G. L. Chiarotti, P. Focher,
!   S. Scandolo, E. Tosatti, and M. Parrinello,
!   J. Phys. Chem. Solids 56, 501 (1995).
!   The function is the one used by ABINIT (input variable ecutsm),
!   see the ABINIT documentation and Appendix B of
!   J. Laflamme Janssen et al., Phys. Rev. B 93, 205147 (2016).
!
!     T = ekin p(x),   x = (emax - ekin) / esmear,
!     p(x) = 1                                  for x > 1,
!     p(x) = 1 / (x^2 (3 + x - 6 x^2 + 3 x^3))   for 0 < x <= 1.
!
!   T and its first and second derivatives are continuous at
!   ekin = emax - esmear, and T diverges at ekin = emax.

subroutine kinetic_smear(ekin, emax, esmear, tkin, dtkin)

! Written 15 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  ekin                            !<  kinetic energy (hartree) of the plane-wave, (1/2)|k+G|^2
  real(REAL64), intent(in)           ::  emax                            !<  kinetic energy cutoff of plane wave expansion (hartree)
  real(REAL64), intent(in)           ::  esmear                          !<  smearing of the cutoff (hartree), no modification if esmear <= 0

! output

  real(REAL64), intent(out)          ::  tkin                            !<  modified kinetic energy (hartree)
  real(REAL64), intent(out)          ::  dtkin                           !<  d tkin / d ekin

! local variables

  real(REAL64)       ::  x                                               !  (emax - ekin) / esmear
  real(REAL64)       ::  d, dd                                           !  denominator 1/p(x) and its derivative

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  XMIN = 1.0E-4_REAL64                      !  avoids the divergence at emax (and beyond, for a fixed basis)


  tkin = ekin
  dtkin = UM

  if(esmear <= ZERO) return

  x = (emax - ekin) / esmear

  if(x >= UM) return

  if(x < XMIN) x = XMIN

  d = x*x*(3 + x*(1 + x*(-6 + 3*x)))
  dd = x*(6 + x*(3 + x*(-24 + 15*x)))

! p = 1/d,  dp/dx = -dd/d^2,  dx/dekin = -1/esmear

  tkin = ekin / d
  dtkin = UM/d + (ekin/esmear) * dd / (d*d)

  return

end subroutine kinetic_smear
