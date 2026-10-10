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

!>  Reduces the metric of a lattice to the metric of a Buerger
!>  (Minkowski in 3D) reduced cell, the cell with the three shortest
!>  non-coplanar lattice vectors.  Real arithmetic version, for use in
!>  geometric searches (see metric_niggli for the exact integer version).
!>
!>  The new lattice vectors are  a'_j = sum_i a_i mtb(i,j),  so that
!>  adotb = mtb^T adot mtb,  det(mtb) = 1, and lattice coordinates
!>  transform as  x = mtb x'.
!>
!>  The reduction sorts the vectors by length, reduces each pair
!>  (b_i <- b_i - nint(b_i.b_j / b_j.b_j) b_j) and replaces b_k by
!>  b_k +- b_i +- b_j when shorter, until nothing changes.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         10 October 2026.
!>  \copyright    GNU Public License v2

subroutine metric_buerger(adot, adotb, mtb)

! Written 10 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

! output

  real(REAL64), intent(out)          ::  adotb(3,3)                      !<  metric of the reduced cell
  integer, intent(out)               ::  mtb(3,3)                        !<  a'_j = sum_i a_i mtb(i,j), det = 1

! local variables

  logical            ::  lchange
  integer            ::  iter
  integer            ::  iq, itmp(3), idet
  real(REAL64)       ::  v(3), vv
  integer            ::  iv(3)

! parameters

  integer, parameter         ::  MXITER = 1000
  real(REAL64), parameter    ::  EPS = 1.0E-12_REAL64

! counters

  integer   ::  i, j, k, si, sj, n


  mtb(:,:) = 0
  do i = 1,3
    mtb(i,i) = 1
  enddo
  adotb(:,:) = adot(:,:)

  do iter = 1,MXITER

    lchange = .FALSE.

!   sorts by length

    do i = 1,2
    do j = i+1,3
      if(adotb(j,j) < adotb(i,i)*(1 - EPS)) then
        itmp(:) = mtb(:,i)
        mtb(:,i) = mtb(:,j)
        mtb(:,j) = itmp(:)
        call metric_buerger_update(adot, mtb, adotb)
      endif
    enddo
    enddo

!   pairwise (Gauss-Lagrange) reduction

    do i = 1,3
    do j = 1,3
      if(i /= j) then
        iq = nint(adotb(i,j) / adotb(j,j))
        if(iq /= 0) then
          mtb(:,i) = mtb(:,i) - iq*mtb(:,j)
          call metric_buerger_update(adot, mtb, adotb)
          lchange = .TRUE.
        endif
      endif
    enddo
    enddo

!   b_k + si b_i + sj b_j

    do k = 1,3
      i = mod(k,3) + 1
      j = mod(k+1,3) + 1
      do si = -1,1,2
      do sj = -1,1,2
        do n = 1,3
          iv(n) = mtb(n,k) + si*mtb(n,i) + sj*mtb(n,j)
        enddo
        do n = 1,3
          v(n) = adot(n,1)*iv(1) + adot(n,2)*iv(2) + adot(n,3)*iv(3)
        enddo
        vv = v(1)*iv(1) + v(2)*iv(2) + v(3)*iv(3)
        if(vv < adotb(k,k)*(1 - EPS)) then
          mtb(:,k) = iv(:)
          call metric_buerger_update(adot, mtb, adotb)
          lchange = .TRUE.
        endif
      enddo
      enddo
    enddo

    if(.not. lchange) exit

  enddo

  if(lchange) then
    write(6,*)
    write(6,'("   STOPPED in metric_buerger:  no convergence after ",    &
         &    i8," iterations")') MXITER
    write(6,*)

    stop

  endif

! makes the determinant positive

  idet = mtb(1,1)*(mtb(2,2)*mtb(3,3) - mtb(2,3)*mtb(3,2))                &
       - mtb(1,2)*(mtb(2,1)*mtb(3,3) - mtb(2,3)*mtb(3,1))                &
       + mtb(1,3)*(mtb(2,1)*mtb(3,2) - mtb(2,2)*mtb(3,1))

  if(idet < 0) then
    mtb(:,:) = -mtb(:,:)
  elseif(idet == 0) then
    write(6,*)
    write(6,'("   STOPPED in metric_buerger:  singular transformation")')
    write(6,*)

    stop

  endif

  return

contains

!>  adotb = mtb^T adot mtb

  subroutine metric_buerger_update(adot, mtb, adotb)

    real(REAL64), intent(in)           ::  adot(3,3)                     !<  metric in direct space
    integer, intent(in)                ::  mtb(3,3)                      !<  a'_j = sum_i a_i mtb(i,j)
    real(REAL64), intent(out)          ::  adotb(3,3)                    !<  metric of the new cell

    integer   ::  i, j, k, l

    do i = 1,3
    do j = 1,3
      adotb(i,j) = 0.0_REAL64
      do k = 1,3
      do l = 1,3
        adotb(i,j) = adotb(i,j) + mtb(k,i)*adot(k,l)*mtb(l,j)
      enddo
      enddo
    enddo
    enddo

    return

  end subroutine metric_buerger_update

end subroutine metric_buerger
