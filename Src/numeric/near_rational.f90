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

!>  Finds if a number is close to a rational or the square-root
!>  of a rational
!>
!>  \author       Jose Luis Martins
!>  \version      4.94
!>  \date         29 May 2014.  13 June 2024.
!>  \copyright    GNU Public License v2

subroutine near_rational(x, lfound, lsquare, nnum, nden, ntry, eps)

! Written 29 May 2014. JLM
! Modified 21 June 2014, bug sqrt(small x). JLM
! Indentation, 6 June 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  real(REAL64), intent(in)           ::  x                               !<  number to be approximated
  integer, intent(in)                ::  ntry                            !<  tries denominators up to ntry
  real(REAL64), intent(in)           ::  eps                             !<  tolerance

! output

  logical, intent(out)               ::  lfound                          !<  an approximation was found
  logical, intent(out)               ::  lsquare                         !<  it is the square of a rational
  integer, intent(out)               ::  nnum                            !<  x ~ nnum/ndem or x**2 ~ nnum/nden and sign of nnum is the same as x
  integer, intent(out)               ::  nden                            !<  x ~ nnum/ndem or x**2 ~ nnum/nden and sign of nnum is the same as x

! local variables

  integer       ::  nmax
  real(REAL64)  ::  xn
  real(REAL64)  ::  eps2

! constants

  real(REAL64), parameter  :: UM = 1.0_REAL64

! counters

  integer       ::  n

  lfound = .FALSE.
  lsquare = .FALSE.
  nnum = 0
  nden = 1

! checks if input makes sense and avoids overflow

  eps2 = max(eps,(UM/2)**30)
  nmax = nint(0.1/eps2)
  nmax = min(max(1,ntry),nmax)

  if(abs(x) < 2.0**30/ntry) then

    do n = 1,nmax

      xn = n*x
      if( abs(xn - nint(xn)) < n*eps2 ) then

        nden = n
        nnum = nint(xn)

        lfound = .TRUE.

        exit

      endif

    enddo

    if(.not. lfound) then

      do n = 1,nmax

        xn = n*x*x
        if(abs(abs(x)-sqrt((nint(xn)*UM)/(n*UM))) < eps2) then

          nden = n
          nnum = nint(xn)

          lfound = .TRUE.
          lsquare = .TRUE.

          exit

        endif

      enddo

    endif

  endif

  return

end subroutine near_rational
