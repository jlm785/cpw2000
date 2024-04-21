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

!>  Calculates traces and eventualy checks symetry of a tensor of the Berry type
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         7 April 2024.
!>  \copyright    GNU Public License v2

subroutine berry_tensor_trace(tensor, trace_n, trace_r, levdegnl, csym, lcheck, mxddeg)

! Written 7 April 2024.  JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  real(REAL64), intent(in)           ::  tensor(mxddeg,mxddeg,3,3)       !<  Berry tensor
  integer, intent(in)                ::  levdegnl                        !<  degeneracy of level
  character(len = 1), intent(in)     ::  csym                            !<  checks symmetry anti-symetry

! output

  real(REAL64), intent(out)          ::  trace_n(3,3)                    !<  Trace on levels of tensor
  real(REAL64), intent(out)          ::  trace_r(mxddeg,mxddeg)          !<  Trace on coordinates of Berry tensor
  logical, intent(out)               ::  lcheck                          !<  result of check

! local variables

  real(REAL64)           ::  xsum


! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-6_REAL64

! counters

  integer    ::  j, k, n


  lcheck = .TRUE.

! trace on levels

  do j = 1,3
  do k = 1,3
    trace_n(j,k) = tensor(1,1,j,k)
    if(levdegnl > 1) then
      do n = 2,levdegnl
        trace_n(j,k) = trace_n(j,k) + tensor(n,n,j,k)
      enddo
    endif
  enddo
  enddo

! trace on coordinates

  do j = 1,levdegnl
  do k = 1,levdegnl
    trace_r(j,k) = tensor(j,k,1,1)
    do n = 2,3
      trace_r(j,k) = trace_r(j,k) + tensor(j,k,n,n)
    enddo
   enddo
  enddo

  if(csym == 'A' .or. csym == 'a') then

    xsum = UM
    do j = 1,3
    do k = 1,3
      xsum = xsum + abs(trace_n(j,k))
    enddo
    enddo
    do j = 1,3
      if(abs(trace_n(j,j)) > xsum*EPS) lcheck = .FALSE.
    enddo
    do j = 1,2
    do k = j+1,3
      if(abs(trace_n(j,k)+trace_n(k,j)) > xsum*EPS) lcheck = .FALSE.
    enddo
    enddo

  elseif(csym == 'S' .or. csym == 's') then

    xsum = UM
    do j = 1,3
    do k = 1,3
      xsum = xsum + abs(trace_n(j,k))
    enddo
    enddo
    do j = 1,2
    do k = j+1,3
      if(abs(trace_n(j,k)-trace_n(k,j)) > xsum*EPS) lcheck = .FALSE.
    enddo
    enddo

  endif


  return

end subroutine berry_tensor_trace
