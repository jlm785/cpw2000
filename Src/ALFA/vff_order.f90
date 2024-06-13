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

!>  Bubble sort by incresing order of the 5 elements of adist
!>  jd and ntau are also permutted with adist.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.11 (1.7 of md)
!>  \date         before January 2013.  5 June 2024.
!>  \copyright    GNU Public License v2

subroutine vff_order(adist, jd, ntau)

! modified January 2013, J.L.Martins, C.S.Loia
! Modified, documentation, details, December 2019. JLM
! Indentation, 5 June 2024. JLM

  implicit none

  integer, parameter      :: REAL64 = selected_real_kind(12)

! input and output

  real(REAL64), intent(inout)   ::  adist(5)                             !<  distance, should be increasing on output
  integer, intent(inout)        ::  jd(5)                                !<  indicates the permutation
  integer, intent(inout)        ::  ntau(3,5)                            !<  permutted tau...

! local variables

  integer            ::  nt
  real(REAL64)       ::  t

! counters

  integer            ::  i, j, k


! Bubble sort is OK for very short arrays

  do i = 1,4
    do j = i+1,5
      if(adist(i) > adist(j)) then
        t = adist(i)
        adist(i) = adist(j)
        adist(j) = t
        nt = jd(i)
        jd(i) = jd(j)
        jd(j) = nt
        do k=1,3
          nt = ntau(k,i)
          ntau(k,i) = ntau(k,j)
          ntau(k,j) = nt
        enddo
      endif
    enddo
  enddo

  return

end subroutine vff_order
