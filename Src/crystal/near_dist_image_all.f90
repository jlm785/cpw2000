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

!>  All the nearest images of a vector taking into account periodicity,
!>  that is all the lattice translations T that minimize |dr + T|,
!>  within a tolerance.  Used to choose the images of the hoppings
!>  in Fourier interpolation (Wigner-Seitz cell).
!>
!>  The caller must provide the Buerger reduced cell (metric_buerger),
!>  so that it is calculated only once for many vectors.  The search is
!>  the one of near_dist_image, and the equidistant images are then
!>  collected in a shell of coefficients -2..2 around the nearest one
!>  in the reduced cell.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         23 October 2026.
!>  \copyright    GNU Public License v2

subroutine near_dist_image_all(nimg, itrans, distmin,                    &
    adotb, mtb, dr, tol,                                                 &
    mxdimg)

! Written 23 October 2026 from near_dist_image. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdimg                          !<  array dimension of the number of images

  real(REAL64), intent(in)           ::  adotb(3,3)                      !<  metric of the Buerger reduced cell
  integer, intent(in)                ::  mtb(3,3)                        !<  a'_j = sum_i a_i mtb(i,j), det = 1, from metric_buerger
  real(REAL64), intent(in)           ::  dr(3)                           !<  vector (lattice coordinates of the original cell)
  real(REAL64), intent(in)           ::  tol                             !<  tolerance for equal distances (same units as the metric, not squared)

! output

  integer, intent(out)               ::  nimg                            !<  number of nearest images
  integer, intent(out)               ::  itrans(3,mxdimg)                !<  dr + itrans(:,n) is the n-th nearest image (lattice coordinates of the original cell)
  real(REAL64), intent(out)          ::  distmin                         !<  distance of the nearest images (not squared)

! local variables

  integer            ::  minv(3,3)                                       !  inverse of mtb (integer, det = 1)
  real(REAL64)       ::  xb(3), zk(3)                                    !  coordinates in the reduced cell
  integer            ::  nfold(3)                                        !  folding of xb
  real(REAL64)       ::  dist, d2min
  integer            ::  kbest(3), ncenter(3), itb(3)
  integer            ::  iter

! parameters

  integer, parameter         ::  MXITER = 100
  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counters

  integer   ::  i, j, k1, k2, k3


! coordinates in the reduced cell, x' = mtb^-1 x

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
      xb(i) = xb(i) + minv(i,j)*dr(j)
    enddo
  enddo
  do i = 1,3
    nfold(i) = nint(xb(i))
    xb(i) = xb(i) - nfold(i)
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
    write(6,'("   STOPPED in near_dist_image_all:  no convergence ",     &
         &    "after ",i8," iterations")') MXITER
    write(6,*)

    stop

  endif

  distmin = sqrt(d2min)

! equidistant images, the translation in the reduced cell is
! ncenter + k - nfold, in the original cell it is mtb times that

  nimg = 0

  do k1 = -2,2
  do k2 = -2,2
  do k3 = -2,2
    zk(1) = xb(1) + ncenter(1) + k1
    zk(2) = xb(2) + ncenter(2) + k2
    zk(3) = xb(3) + ncenter(3) + k3
    dist = ZERO
    do i = 1,3
    do j = 1,3
      dist = dist + zk(i)*adotb(i,j)*zk(j)
    enddo
    enddo
    dist = sqrt(max(dist,ZERO))

    if(dist - distmin < tol) then

      nimg = nimg + 1

      if(nimg > mxdimg) then
        write(6,*)
        write(6,'("   STOPPED in near_dist_image_all:  more than ",      &
             &    i5," images, increase mxdimg")') mxdimg
        write(6,*)

        stop

      endif

      itb(1) = ncenter(1) + k1 - nfold(1)
      itb(2) = ncenter(2) + k2 - nfold(2)
      itb(3) = ncenter(3) + k3 - nfold(3)
      do i = 1,3
        itrans(i,nimg) = mtb(i,1)*itb(1) + mtb(i,2)*itb(2) + mtb(i,3)*itb(3)
      enddo

    endif

  enddo
  enddo
  enddo

  return

end subroutine near_dist_image_all
