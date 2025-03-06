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

!>  Performs a shift of the vectors that were indicated as linearly
!>  dependent by irow.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         August 17 2015, March 5 2025.
!>  \copyright    GNU Public License v2

subroutine grsch_shift_c16(xvec, mtxd, nconv, nvec, irow, ired,          &
                           mxddim)

! Written August 17 2015. JLM
! Modified, documentation, January 2020. JLM
! Modified, indentation, March 5 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  nconv                           !<  number of vectors that are already orthogonal
  integer, intent(in)                ::  nvec                            !<  number of vectors

  integer, intent(in)                ::  irow(nvec)                      !<  if irow(i)=0 the vector is linearly dependent on the others and the corresponding xvec is zero

! input and output

  complex(REAL64), intent(inout)     ::  xvec(mxddim,nvec)               !<  component j of vector i

! output

  integer, intent(out)               ::  ired                            !<  number of vectors eliminated (0s of irow)

! counters

  integer    ::   n, j

  j = nconv+1
  ired = 0

  do n = nconv+1,nvec
    if(irow(n) == 0) then
      ired = ired + 1
    else
      if(j < n) then
        call zcopy(mtxd, xvec(1,n), 1, xvec(1,j), 1)
      endif
      j = j + 1
    endif
  enddo

  return

end subroutine grsch_shift_c16

