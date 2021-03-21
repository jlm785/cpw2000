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

!>  Calculates the product np for semiconductors multiplied by exp(E_gap/k_B T)
!>  Assumes chemical potential is not near the band edges

subroutine opt_rad_np(xnp, egap, tau, adot, evbmez, ehist, dhist, nhist)

! Written 16 October 2020. JLM
! copyright  Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter                  ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                 ::  nhist                          !<  number of points in histogram

  real(REAL64), intent(in)            ::  tau                            !<  temperature in au
  
  real(REAL64), intent(in)            ::  adot(3,3)                      !<  metric in real space
  real(REAL64), intent(in)            ::  evbmez                         !<  cond. band max should be above this value
  real(REAL64), intent(in)            ::  ehist(nhist)                   !<  histogram energy
  real(REAL64), intent(in)            ::  dhist(nhist)                   !<  histogram dos

! output

  real(REAL64), intent(out)           ::  xnp                            !<  n*p * exp(E_gap/k_B T)
  real(REAL64), intent(out)           ::  egap                           !<  band gap (underestimate)

! local variables

  real(REAL64)       ::  evbm                                            !  energy of maximum of valence band
  real(REAL64)       ::  ecbm                                            !  energy of minimum of conduction band

  real(REAL64)       ::  xp                                              !  density of holes without exp(E-mu), atomic units
  real(REAL64)       ::  xn                                              !  density of electrons without exp(mu-E), atomic units

  integer            ::  icbm, ivbm                                      !  indices associated with CBM and VBM
  integer            ::  istart

  real(REAL64)       ::  vcell, bdot(3,3)

  logical            ::  lmetal

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64
  real(REAL64), parameter ::  EPS = 1.0E-10_REAL64

! counters

  integer            ::  i

! identifies band edges according to DOS, not the "rough" estimates

  do i = 1,nhist
    if(ehist(i) > evbmez) then
      istart = i

      exit

    endif
  enddo

  if(abs(dhist(i)) > EPS) then

    write(6,'("   error in opt_rad_np :  mid of the gap misidentified.")')

    egap = ZERO
    xnp = ZERO

    lmetal = .TRUE.

  endif

  if(.not. lmetal) then
    do i = istart, 1,-1
      ivbm = i
      evbm = ehist(i)

      if(abs(dhist(i)) > EPS) exit

    enddo

    do i = istart, nhist
      icbm = i
      ecbm = ehist(i)

      if(abs(dhist(i)) > EPS) exit

    enddo

    if(ivbm > icbm) then
      write(6,'("  opt_rad_np:  ivbm, icbm = ",2i8)')

      stop

    endif

    egap = ecbm - evbm

! ! fast integration, but not too accurate

    xp = ZERO

    do i = 1,ivbm-1
      xp = xp + (ehist(i+1)-ehist(i))*(dhist(i+1)*exp((ehist(i+1)-evbm)/tau)+dhist(i)*exp((ehist(i)-evbm)/tau))/2
    enddo

    xn = ZERO

    do i = icbm,nhist-1
      xn = xn + (ehist(i+1)-ehist(i))*(dhist(i+1)*exp((ecbm-ehist(i+1))/tau)+dhist(i)*exp((ecbm-ehist(i))/tau))/2
    enddo

    call adot_to_bdot(adot,vcell,bdot)

    xn = xn/vcell
    xp = xp/vcell

    xnp = xn*xp

  endif


  return
end subroutine opt_rad_np
