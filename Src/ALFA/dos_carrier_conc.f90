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

!>  Finds the carrier concentration for a semiconductor
!>  also gives their derivatives with respect to chemical potential

subroutine dos_carrier_conc(xn, xp, dxndef, dxpdef, tau, ef,              &
                            nvbm, ncbm, nhist, ehist, dhist)

! written 5 December 2013.
! copyright Jose Luis Martins / INESC-MN
! Modified, documentation, 19 September 2020. JLM
! Modified, limits of integration, 19 October 2020. JLM

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)         ::  tau                                !<  temperature in atomic units
  real(REAL64), intent(in)         ::  ef                                 !<  Fermi energy in atomic units
  integer, intent(in)              ::  nvbm                               !<  ehist(nvbm) is above valence band maximum
  integer, intent(in)              ::  ncbm                               !<  ehist(ncbm) is below conduction band minimum
  integer, intent(in)              ::  nhist                              !<  size of energy array
  real(REAL64), intent(in)         ::  ehist(nhist)                       !<  energy array (Hartree)
  real(REAL64), intent(in)         ::  dhist(nhist)                       !<  density of states array (electron/cell/Hartree)

! output

  real(REAL64), intent(out)        ::  xn                                 !<  concentration of electrons (per unit cell)
  real(REAL64), intent(out)        ::  xp                                 !<  concentration of holes (per unit cell)
  real(REAL64), intent(out)        ::  dxndef                             !<  d xn / d ef
  real(REAL64), intent(out)        ::  dxpdef                             !<  d xp / d ef

! local variables

  real(REAL64)    :: arg, ff, dffdef, argref

! counters

  integer               ::  n

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! use trapezoidal rule as the data maybe noisy

  argref = (ef - ehist(nvbm)) / tau
  xp = ZERO
  dxpdef = ZERO
  do n = nvbm,2,-1
    arg = (ef - ehist(n)) / tau
    call fermi_dirac(arg,ff,dffdef)
    dffdef = dffdef / tau
    if(arg - argref > 35 ) then

      exit

    endif
    xp = xp + ff*dhist(n)*(ehist(n+1)-ehist(n-1)) / 2
    dxpdef = dxpdef + dffdef*dhist(n)*(ehist(n+1)-ehist(n-1)) / 2
  enddo

  argref = (ehist(ncbm) - ef) / tau
  xn = ZERO
  dxndef = ZERO
  do n = ncbm,nhist-1
    arg = (ehist(n) - ef) / tau
    call fermi_dirac(arg,ff,dffdef)
    dffdef = -dffdef / tau
    if(arg - argref > 35 ) then

      exit

    endif
    xn = xn + ff*dhist(n)*(ehist(n+1)-ehist(n-1)) / 2
    dxndef = dxndef + dffdef*dhist(n)*(ehist(n+1)-ehist(n-1)) / 2
  enddo

  return
end subroutine dos_carrier_conc
