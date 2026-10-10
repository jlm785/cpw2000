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

!>  Closest distance between two atoms taking into account periodicity,
!>  and the corresponding image of r1 - r0 (nearest image vector).
!>
!>  The search is done in the Buerger (Minkowski) reduced cell.  The
!>  difference of coordinates is folded with nint and the 26 neighbours
!>  with coefficients -1, 0, 1 are searched, recentering on the best one
!>  until the center is the best.  For a reduced cell in 3D all
!>  Voronoi-relevant vectors have coefficients -1, 0, 1, so the result
!>  is the exact nearest image for any cell, however oblique.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         10 October 2026.
!>  \copyright    GNU Public License v2

subroutine near_dist_image(distmin, rmin, adot, r0, r1)

! Written 10 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  r0(3)                           !<  position of the first atom (lattice coordinates)
  real(REAL64), intent(in)           ::  r1(3)                           !<  position of the second atom (lattice coordinates)

! output

  real(REAL64), intent(out)          ::  distmin                         !<  closest distance between images (not squared)
  real(REAL64), intent(out)          ::  rmin(3)                         !<  nearest image of r1 - r0 (lattice coordinates)

! local variables

  real(REAL64)       ::  adotb(3,3)                                      !  metric of the reduced cell
  integer            ::  mtb(3,3)                                        !  a'_j = sum_i a_i mtb(i,j)
  integer            ::  minv(3,3)                                       !  inverse of mtb (integer, det = 1)
  real(REAL64)       ::  xb(3), zk(3)                                    !  coordinates in the reduced cell
  real(REAL64)       ::  dist, d2min
  integer            ::  kbest(3), ncenter(3)
  integer            ::  iter

! parameters

  integer, parameter         ::  MXITER = 100
  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counters

  integer   ::  i, j, k1, k2, k3


! Buerger reduced cell and coordinates in that cell, x' = mtb^-1 x

  call metric_buerger(adot, adotb, mtb)

  minv(1,1) = mtb(2,2)*mtb(3,3) - mtb(2,3)*mtb(3,2)
  minv(1,2) = mtb(1,3)*mtb(3,2) - mtb(1,2)*mtb(3,3)
  minv(1,3) = mtb(1,2)*mtb(2,3) - mtb(1,3)*mtb(2,2)
  minv(2,1) = mtb(2,3)*mtb(3,1) - mtb(2,1)*mtb(3,3)
  minv(2,2) = mtb(1,1)*mtb(3,3) - mtb(1,3)*mtb(3,1)
  minv(2,3) = mtb(1,3)*mtb(2,1) - mtb(1,1)*mtb(2,3)
  minv(3,1) = mtb(2,1)*mtb(3,2) - mtb(2,2)*mtb(3,1)
  minv(3,2) = mtb(1,2)*mtb(3,1) - mtb(1,1)*mtb(3,2)
  minv(3,3) = mtb(1,1)*mtb(2,2) - mtb(1,2)*mtb(2,1)

  do i = 1,3
    xb(i) = ZERO
    do j = 1,3
      xb(i) = xb(i) + minv(i,j)*(r1(j) - r0(j))
    enddo
  enddo
  do i = 1,3
    xb(i) = xb(i) - nint(xb(i))
  enddo

! search of the 26 neighbours, recentering until the center is the best

  ncenter(:) = 0
  d2min = ZERO
  do i = 1,3
  do j = 1,3
    d2min = d2min + xb(i)*adotb(i,j)*xb(j)
  enddo
  enddo

  do iter = 1,MXITER

    kbest(:) = 0

    do k1 = -1,1
    do k2 = -1,1
    do k3 = -1,1
      zk(1) = xb(1) + ncenter(1) + k1
      zk(2) = xb(2) + ncenter(2) + k2
      zk(3) = xb(3) + ncenter(3) + k3
      dist = ZERO
      do i = 1,3
      do j = 1,3
        dist = dist + zk(i)*adotb(i,j)*zk(j)
      enddo
      enddo
      if(dist < d2min) then
        d2min = dist
        kbest(1) = k1
        kbest(2) = k2
        kbest(3) = k3
      endif
    enddo
    enddo
    enddo

    if(kbest(1) == 0 .and. kbest(2) == 0 .and. kbest(3) == 0) exit

    ncenter(:) = ncenter(:) + kbest(:)

  enddo

  if(iter > MXITER) then
    write(6,*)
    write(6,'("   STOPPED in near_dist_image:  no convergence after ",   &
         &    i8," iterations")') MXITER
    write(6,*)

    stop

  endif

  distmin = sqrt(d2min)

! back to the original lattice coordinates, x = mtb x'

  do i = 1,3
    zk(i) = xb(i) + ncenter(i)
  enddo
  do i = 1,3
    rmin(i) = ZERO
    do j = 1,3
      rmin(i) = rmin(i) + mtb(i,j)*zk(j)
    enddo
  enddo

  return

end subroutine near_dist_image
