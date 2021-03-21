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

!>  finds a reasonable scale for a density optical response plot

  subroutine opt_plot_range(el, neig, nval, nhtarg, deltae, nhist, nrk, mxdbnd)

! Adapted from dos_plot_range, 11 December 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands
  integer, intent(in)                ::  nrk                             !<  number of irreducible k-points

  real(REAL64), intent(in)           ::  el(mxdbnd,nrk)                  !<  eigenvalues in Hartree

  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  integer, intent(in)                ::  nval                            !<  number of valence bands

  integer, intent(in)                ::  nhtarg                          !<  approximate number of points for plot

! output

  real(REAL64), intent(out)          ::  deltae                          !<  recomended energy step for histogram
  integer, intent(out)               ::  nhist                           !<  recomended number of points for histogram

! local variables

  real(REAL64)                 ::  emax, b

! constants

  real(REAL64), parameter      ::  HARTREE = 27.21138386_REAL64
  real(REAL64), parameter      ::  ZERO = 0.0_REAL64

! counters

  integer   ::  irk


! find emax Maximum rande, minimum difference)

  emax = el(neig,1) - el(nval,1)
  do irk = 1,nrk
    if(emax < el(neig,irk) - el(nval,irk)) emax = el(neig,irk) - el(nval,irk)
  enddo

! find scale for the plot

  deltae = HARTREE*emax / nhtarg

  call plot_step(deltae,b)

  deltae = b/HARTREE

  emax = emax + deltae
  nhist = nint(emax/deltae) + 1

  return
end subroutine opt_plot_range
