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

!>  Resets the number of states according to degeneracy.
!>  neigin <= neig <= neigtmp
!>  ei(neig) and ei(neig+1) should not be degenerate.
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         24 October 2023.
!>  \copyright    GNU Public License v2

subroutine berry_reset_neig(ipr, neig, neigin, neigtmp, ei,              &
             mxdbnd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  ipr                             !<  if not zero prints comments on the result

  integer, intent(in)                ::  neigin                          !<  minimum desired number of bands
  integer, intent(in)                ::  neigtmp                         !<  number of bands for which energies were calculated

  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue in hartree

! output

  integer, intent(out)               ::  neig                            !<  number of bands taking into account degeneracy

! local variables

  real(REAL64)              ::  tol

! parameters

  real(REAL64), parameter       ::  EPS = 1.0E-10_REAL64

! counters

  integer    ::  n


! paranoid checks

  if(neigin > mxdbnd .or. neigtmp > mxdbnd .or. neigin > neigtmp) then
    write(6,*)
    write(6,*) '   STOPPED in berry_reset_neig.  Inconsisten input.'
    write(6,*) '   neigin = ', neigin, '  neigtmp = ', neigtmp, '  mxdbnd = ', mxdbnd
    write(6,*)

    stop

  endif

  if(neigtmp > 1) then
    do n = 2,neigtmp
      if(ei(n) + EPS < ei(n-1)) then
        write(6,*)
        write(6,*) '   STOPPED in berry_reset_neig.  eigenvalues are not increasing.'
        write(6,*) '   n = ', n, '  ei(n) = ', ei(n), '  ei(n-1) = ', ei(n-1)
        write(6,*)

        stop

      endif
    enddo
  endif

! identifies desired number of eigenvectors

  if(neigin == neigtmp) then
    neig = neigin
  else
    tol = (ei(neigtmp) - ei(1))/1000 +10*EPS
    do n = neigin, neigtmp-1
      neig = n

      if(ei(n+1) - ei(n) > tol) exit

    enddo
  endif

! prints eventual comments

  if(ipr /= 0) then
    write(6,*)

    if(neigin == neigtmp) then
      write(6,*) '  WARNING in berry_reset_neig. Not enough bands'
      write(6,*) '  neigin = neigtmp = neig = ', neig
    endif
    if(neig == neigtmp) then
      write(6,*) '  WARNING in berry_reset_neig. Not enough bands'
      write(6,*) '  neigin = ', neigin,'  neigtmp = neig = ', neig
    endif
    if(neig > neigin) then
      write(6,*) '  Number of eigenvalues increased from ', neigin, ' to ',neig
    endif

    write(6,*)
  endif

  return

end subroutine berry_reset_neig




