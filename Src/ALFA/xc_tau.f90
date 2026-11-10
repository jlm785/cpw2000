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

!>  Computes the Thomas-Fermi von Weizsaker kinetic energy
!   functional and its derivatives.
!>
!>  \author       J.L. Martins
!>  \version      5.13
!>  \date         22 April 2026.
!>  \copyright    GNU Public License v2

subroutine xc_tau(rho, grho, tauunif, tausingle,                         &
                  d_tauunif_dr, d_tausingle_dr, d_tausingle_dgr)

! Written 22 April 2026, extracted xc_mgga. JLM
! Moved from xc_mgga.f90 to its own file, 10 November 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)

! output

  real(REAL64), intent(out)          ::  tauunif                         !<  Thomas-Fermi kinetic energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  tausingle                       !<  von Weizsaker correction to tauunif (hartree/bohr^3)

  real(REAL64), intent(out)          ::  d_tauunif_dr                    !<  d tauunif / d rho
  real(REAL64), intent(out)          ::  d_tausingle_dr                  !<  d tausingle / d rho
  real(REAL64), intent(out)          ::  d_tausingle_dgr                 !<  d tausingle / d grho

! local variables

  real(REAL64)     ::  rs                                                !  Wigner-Seitz radius
  real(REAL64)     ::  xkf, d_xkf_dr                                     !  Fermi wave-vector and derivative
  real(REAL64)     ::  grloc                                             !  local value of grho

! parameters

  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  AKF = (9*PI/4)**(UM/(3*UM))


  rs = (3/(4*PI*rho))**(UM/3)

! avoids noisy data

  grloc = grho
  if(grho < ZERO) grloc = ZERO


  xkf = AKF / rs                                                         !  (3*PI*PI * rho)**(UM/3)
  d_xkf_dr = (UM/3) * xkf / rho                                          !  (UM/3)*(3*PI*PI * rho)**(UM/3) / rho

  tauunif = (3*UM / 10) * xkf*xkf * rho
  d_tauunif_dr = (3*UM / 10) * xkf*xkf + 2*(3*UM / 10) * xkf*d_xkf_dr * rho

  tausingle = grloc*grloc / (8*rho)                                      !   xkf*xkf * rho * s*s /2
  d_tausingle_dr = -tausingle/rho
  d_tausingle_dgr = 2* grloc / (8*rho)

  return

end subroutine xc_tau
