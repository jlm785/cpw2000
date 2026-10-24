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

!>  Weights of the squared errors |H psi - e psi|^2 of the eigenvectors,
!>  so that the occupied states are converged to epspsi and the
!>  empty states to a looser tolerance:
!>
!   states 1 to nocc                         epspsi
!   state nocc+1                             3 epspsi
!   states nocc+2 to nocc + (neig-nocc)/2    10 epspsi
!   higher states                            100 epspsi
!
!   The loose tolerances are capped at max(epspsi,EPSCAP).  The weight
!   is (epspsi/tolerance)**2, and the weighted errors are compared
!   with epspsi**2.  If nocc >= neig all states (including those added
!   for degeneracies beyond neig) have weight one.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         24 October 2026.
!>  \copyright    GNU Public License v2

subroutine ditsp_weight(nocc, neig, epspsi, wgt, mxdbnd)

! Written 24 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  nocc                            !<  number of occupied states (converged to epspsi)
  integer, intent(in)                ::  neig                            !<  number of requested eigenvectors
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the occupied eigenvectors

! output

  real(REAL64), intent(out)          ::  wgt(mxdbnd)                     !<  weight of the squared error of eigenvector n

! local variables

  integer           ::  nocce                                            !  nocc between 1 and neig
  real(REAL64)      ::  tolmax                                           !  loosest tolerance
  real(REAL64)      ::  tol

! parameters

  real(REAL64), parameter     ::  UM = 1.0_REAL64
  real(REAL64), parameter     ::  EPSCAP = 1.0E-3_REAL64                 !  cap of the loose tolerances

! counters

  integer   ::  n


  nocce = max(1,min(nocc,neig))

  if(nocce >= neig) then

    do n = 1,mxdbnd
      wgt(n) = UM
    enddo

  else

    tolmax = max(epspsi,EPSCAP)

    do n = 1,mxdbnd
      if(n <= nocce) then
        tol = epspsi
      elseif(n == nocce+1) then
        tol = 3*epspsi
      elseif(n <= nocce + (neig-nocce)/2) then
        tol = 10*epspsi
      else
        tol = 100*epspsi
      endif
      tol = min(tol,tolmax)
      wgt(n) = (epspsi/tol)*(epspsi/tol)
    enddo

  endif

  return

end subroutine ditsp_weight
