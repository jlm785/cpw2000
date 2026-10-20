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

!>  Solves the assignment problem: finds the one-to-one
!>  assignment of rows to columns that maximizes the sum of
!>  score(i,iassign(i)).  Hungarian algorithm with shortest
!>  augmenting paths, O(n**3) with n = max(nrow,ncol).
!>  The matrix is padded with zeros if it is not square,
!>  so if nrow > ncol some rows are left unassigned.
!>
!   References:
!   H. W. Kuhn, Naval Research Logistics Quarterly 2, 83 (1955).
!   J. Munkres, J. SIAM 5, 32 (1957).
!   R. Jonker and A. Volgenant, Computing 38, 325 (1987).
!   The structure (potentials u, v and the augmenting path way)
!   follows the version in https://cp-algorithms.com,
!   "Hungarian algorithm for solving the assignment problem".
!>
!>  \author       claude
!>  \version      5.13
!>  \date         19 October 2026.
!>  \copyright    GNU Public License v2

subroutine num_assign_max(nrow, ncol, score, iassign, mxdrow)

! Written 19 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdrow                          !<  array dimension of rows

  integer, intent(in)                ::  nrow                            !<  number of rows
  integer, intent(in)                ::  ncol                            !<  number of columns
  real(REAL64), intent(in)           ::  score(mxdrow,ncol)              !<  score of assigning row i to column j

! output

  integer, intent(out)               ::  iassign(nrow)                   !<  column assigned to row i, 0 if none

! local allocatable arrays

  real(REAL64), allocatable          ::  cost(:,:)                       !  square cost matrix (minus the score)
  real(REAL64), allocatable          ::  u(:), v(:)                      !  dual variables (potentials) of rows and columns
  real(REAL64), allocatable          ::  xminv(:)                        !  minimum reduced cost reaching column j
  integer, allocatable               ::  irowcol(:)                      !  row assigned to column j (0 is the auxiliary column)
  integer, allocatable               ::  iway(:)                         !  previous column in the augmenting path
  logical, allocatable               ::  lused(:)                        !  column already in the path

! local variables

  integer           ::  n
  integer           ::  j0, j1, i0
  real(REAL64)      ::  delta, cur

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  real(REAL64), parameter     ::  XBIG = 1.0E+30_REAL64

! counters

  integer   ::  i, j


  n = max(nrow,ncol)

  do i = 1,nrow
    iassign(i) = 0
  enddo

  if(n < 1) return

  allocate(cost(n,n))
  allocate(u(0:n), v(0:n), xminv(0:n))
  allocate(irowcol(0:n), iway(0:n))
  allocate(lused(0:n))

  do j = 1,n
  do i = 1,n
    cost(i,j) = ZERO
    if(i <= nrow .and. j <= ncol) cost(i,j) = -score(i,j)
  enddo
  enddo

  do j = 0,n
    u(j) = ZERO
    v(j) = ZERO
    irowcol(j) = 0
    iway(j) = 0
  enddo

! adds the rows one by one, each with an augmenting path

  do i = 1,n

    irowcol(0) = i
    j0 = 0
    do j = 0,n
      xminv(j) = XBIG
      lused(j) = .FALSE.
    enddo

    do
      lused(j0) = .TRUE.
      i0 = irowcol(j0)
      delta = XBIG
      j1 = 0
      do j = 1,n
        if(.not. lused(j)) then
          cur = cost(i0,j) - u(i0) - v(j)
          if(cur < xminv(j)) then
            xminv(j) = cur
            iway(j) = j0
          endif
          if(xminv(j) < delta) then
            delta = xminv(j)
            j1 = j
          endif
        endif
      enddo

      do j = 0,n
        if(lused(j)) then
          u(irowcol(j)) = u(irowcol(j)) + delta
          v(j) = v(j) - delta
        else
          xminv(j) = xminv(j) - delta
        endif
      enddo

      j0 = j1
      if(irowcol(j0) == 0) exit

    enddo

!   inverts the augmenting path

    do
      j1 = iway(j0)
      irowcol(j0) = irowcol(j1)
      j0 = j1
      if(j0 == 0) exit
    enddo

  enddo

  do j = 1,ncol
    if(irowcol(j) >= 1 .and. irowcol(j) <= nrow) iassign(irowcol(j)) = j
  enddo

  deallocate(cost)
  deallocate(u, v, xminv)
  deallocate(irowcol, iway)
  deallocate(lused)

  return

end subroutine num_assign_max
