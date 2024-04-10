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

!>  Calculates and eventualy checks symetry of a tensor of the Berry type
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         7 April 2024.
!>  \copyright    GNU Public License v2

subroutine berry_tensor_trace(tensor, ttrace, pvec, levdegnl, check, lcheck, mxddeg)

! Written 7 April 2024.  JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  real(REAL64), intent(in)           ::  tensor(mxddeg,mxddeg,3,3)       !<  Berry tensor
  integer, intent(in)                ::  levdegnl                        !<  degeneracy of level
  character(len = 1), intent(in)     ::  check                           !<  checks symmetry anti-symetry

! output

  real(REAL64), intent(out)          ::  ttrace(3,3)                     !<  Trace of Berry tensor
  real(REAL64), intent(out)          ::  pvec(3)                         !<  if check = 'A' calculates the associated pseudo-vector.
  logical, intent(out)               ::  lcheck                          !<  result of check

! local variables

  real(REAL64)           ::  xsum


! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-6_REAL64

! counters

  integer    ::  j, k, n


  pvec(:) = ZERO
  lcheck = .TRUE.

  do j = 1,3
  do k = 1,3
    ttrace(j,k) = tensor(1,1,j,k)
    if(levdegnl > 1) then
      do n = 2,levdegnl
        ttrace(j,k) = ttrace(j,k) + tensor(n,n,j,k)
      enddo
    endif
  enddo
  enddo

  if(check == 'A' .or. check == 'a') then

    xsum = UM
    do j = 1,3
    do k = 1,3
      xsum = xsum + abs(ttrace(j,k))
    enddo
    enddo
    do j = 1,3
      if(abs(ttrace(j,j)) > xsum*EPS) lcheck = .FALSE.
    enddo
    do j = 1,2
    do k = j+1,3
      if(abs(ttrace(j,k)+ttrace(k,j)) > xsum*EPS) lcheck = .FALSE.
    enddo
    enddo

    if(lcheck) then
      pvec(1) = ttrace(2,3)
      pvec(2) = ttrace(3,1)
      pvec(3) = ttrace(1,2)
    endif

  elseif(check == 'S' .or. check == 's') then

    xsum = UM
    do j = 1,3
    do k = 1,3
      xsum = xsum + abs(ttrace(j,k))
    enddo
    enddo
    do j = 1,2
    do k = j+1,3
      if(abs(ttrace(j,k)-ttrace(k,j)) > xsum*EPS) lcheck = .FALSE.
    enddo
    enddo

  endif


  return

end subroutine berry_tensor_trace
