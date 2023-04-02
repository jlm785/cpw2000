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

!>  Calculates the auto-correlation function
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         5 March 2023.
!>  \copyright    GNU Public License v2

subroutine plot_z1D_auto_corr(nn, ave, kstart, fac, xpeak, lfound)

! extracted from plot_rho_v_average, 2 March 2023. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nn                              !<  number of points

  real(REAL64), intent(in)           ::  ave(nn)                         !<  layer average of function in the 3d direction of fft grid
  integer, intent(in)                ::  kstart                          !<  starting point to search for a maximum

! output

  real(REAL64), intent(out)          ::  fac(0:nn/2)                     !<  autocorrelation function
  real(REAL64), intent(out)          ::  xpeak                           !<  position of nearest maximum
  logical, intent(out)               ::  lfound                          !<  true if a peak was found away from extrema

! local variables

  integer              ::  ku, kd, ipeak
  real(REAL64)         ::  secder

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer       ::  i, j, k



!   quick and dirty auto-correlation function. See Wiener–Khinchin theorem

    do k = 0,nn/2
      fac(k) = ZERO
      do j = 1,nn
        i = mod(j + k-1,nn) + 1
        fac(k) = fac(k) + ave(j)*ave(i)
      enddo
      fac(k) = fac(k)/nn
    enddo

!   searches from the average layer distance

    ipeak = 0
    do k = 1, nn/2
      ku = kstart + k - 1
      kd = kstart - k

      if(0 < ku .and. ku < nn/2) then
        if(fac(ku) > fac(ku-1) .and. fac(ku) > fac(ku+1) ) then
          ipeak = ku

          exit

        endif
      endif

      if(0 < kd .and. kd < nn/2) then
        if(fac(kd) > fac(kd-1) .and. fac(kd) > fac(kd+1) ) then
          ipeak = kd

          exit

        endif
      endif
    enddo

!   second order interpolation around the maximum

    lfound = .FALSE.
    if(ipeak < 1) then
      xpeak = ZERO
    elseif(ipeak > nn/2-1) then
      xpeak = UM*(nn/2)
    else
      secder = fac(ipeak-1) - 2*fac(ipeak) + fac(ipeak+1)
      if(secder < ZERO) then
        xpeak = ipeak*UM - (fac(ipeak+1) - fac(ipeak-1)) / (2*secder)
        if(xpeak > (ipeak-1)*UM .and. xpeak < (ipeak+1)*UM) lfound = .TRUE.
      else
        xpeak = ipeak*UM
      endif
    endif

  return

end subroutine plot_z1D_auto_corr
