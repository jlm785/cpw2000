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

!>  Converts a general Monkhorst-Pack mesh (block kgrid_Monkhorst_Pack,
!>  with off-diagonal elements as in siesta) into a box mesh nx*ny*nz
!>  with shifts sx,sy,sz and an integer unimodular matrix kmat.
!>
!>
!>  \author       claude
!>  \version      5.13
!>  \date         19 November 2026.
!>  \copyright    GNU Public License v2

subroutine kgrid_box(mpgrid, mpshift, nx, ny, nz, sx, sy, sz, kmat)

! Written 19 November 2026. claude

! The supercell vectors are  A_i = sum_j mpgrid(i,j) a_j  and the mesh
! points are  k = sum_i (n_i + mpshift(i)) B_i,  B_i the reciprocal
! vectors of the supercell.  In reciprocal lattice coordinates
! rk = mpgrid^-1 (n + mpshift).  With integer row and column operations
! L mpgrid R = D = diag(nx,ny,nz)  (diagonal, not necessarily the Smith
! normal form), and
!
! rk = kmat ( ((i-1)+sx)/nx, ((j-1)+sy)/ny, ((k-1)+sz)/nz )
!
! with  kmat = R  and  (sx,sy,sz) = L mpshift  (modulo 1).
! A diagonal mpgrid with positive elements gives kmat = identity.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mpgrid(3,3)                     !<  supercell of the mesh,  A_i = sum_j mpgrid(i,j) a_j  (line i of kgrid_Monkhorst_Pack)
  real(REAL64), intent(in)           ::  mpshift(3)                      !<  shift of the mesh in units of the reciprocal vectors of the supercell (usually 0.5 or 0.0)

! output

  integer, intent(out)               ::  nx                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 1
  integer, intent(out)               ::  ny                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 2
  integer, intent(out)               ::  nz                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 3
  real(REAL64), intent(out)          ::  sx                              !<  offset of the integration mesh, direction 1
  real(REAL64), intent(out)          ::  sy                              !<  offset of the integration mesh, direction 2
  real(REAL64), intent(out)          ::  sz                              !<  offset of the integration mesh, direction 3
  integer, intent(out)               ::  kmat(3,3)                       !<  rk = kmat (((i-1)+sx)/nx, ((j-1)+sy)/ny, ((k-1)+sz)/nz), det = +-1

! local variables

  integer          ::  a(3,3)                                            !  working copy, L mpgrid R
  integer          ::  lmat(3,3)                                         !  row operations
  integer          ::  ip, jp, amin, iq, itmp
  logical          ::  ldiag, lclean
  real(REAL64)     ::  s(3)

! parameters

  integer, parameter         ::  MXITER = 1000
  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-10_REAL64

! counters

  integer   ::  i, j, it, iter


  ldiag = .TRUE.
  do j = 1,3
  do i = 1,3
    if(i /= j .and. mpgrid(i,j) /= 0) ldiag = .FALSE.
  enddo
  if(mpgrid(j,j) < 1) ldiag = .FALSE.
  enddo

  do j = 1,3
  do i = 1,3
    kmat(i,j) = 0
    lmat(i,j) = 0
    a(i,j) = mpgrid(i,j)
  enddo
  kmat(j,j) = 1
  lmat(j,j) = 1
  enddo

  if(.not. ldiag) then

    do it = 1,3

      do iter = 1,MXITER

!       pivot, the smallest non-zero element of the remaining block (diagonal first)

        ip = 0
        jp = 0
        amin = 0
        if(a(it,it) /= 0) then
          ip = it
          jp = it
          amin = abs(a(it,it))
        endif
        do j = it,3
        do i = it,3
          if(a(i,j) /= 0 .and. (amin == 0 .or. abs(a(i,j)) < amin)) then
            ip = i
            jp = j
            amin = abs(a(i,j))
          endif
        enddo
        enddo

        if(amin == 0) then
          write(6,*)
          write(6,*) '   STOPPED in kgrid_box:  the matrix of the',      &
                     ' Monkhorst-Pack mesh is singular'
          write(6,'(3x,3i6)') (mpgrid(1,j),j=1,3)
          write(6,'(3x,3i6)') (mpgrid(2,j),j=1,3)
          write(6,'(3x,3i6)') (mpgrid(3,j),j=1,3)
          write(6,*)

          stop

        endif

!       moves the pivot to (it,it)

        if(ip /= it) then
          do j = 1,3
            itmp = a(it,j)
            a(it,j) = a(ip,j)
            a(ip,j) = itmp
            itmp = lmat(it,j)
            lmat(it,j) = lmat(ip,j)
            lmat(ip,j) = itmp
          enddo
        endif
        if(jp /= it) then
          do i = 1,3
            itmp = a(i,it)
            a(i,it) = a(i,jp)
            a(i,jp) = itmp
            itmp = kmat(i,it)
            kmat(i,it) = kmat(i,jp)
            kmat(i,jp) = itmp
          enddo
        endif

!       reduces the column (row operations) and the row (column operations)

        do i = it+1,3
          iq = a(i,it) / a(it,it)
          if(iq /= 0) then
            do j = 1,3
              a(i,j) = a(i,j) - iq*a(it,j)
              lmat(i,j) = lmat(i,j) - iq*lmat(it,j)
            enddo
          endif
        enddo
        do j = it+1,3
          iq = a(it,j) / a(it,it)
          if(iq /= 0) then
            do i = 1,3
              a(i,j) = a(i,j) - iq*a(i,it)
              kmat(i,j) = kmat(i,j) - iq*kmat(i,it)
            enddo
          endif
        enddo

        lclean = .TRUE.
        do i = it+1,3
          if(a(i,it) /= 0 .or. a(it,i) /= 0) lclean = .FALSE.
        enddo

        if(lclean) exit

      enddo

      if(.not. lclean) then
        write(6,*)
        write(6,*) '   STOPPED in kgrid_box:  no convergence'
        write(6,*)

        stop

      endif

      if(a(it,it) < 0) then
        do j = 1,3
          a(it,j) = -a(it,j)
          lmat(it,j) = -lmat(it,j)
        enddo
      endif

    enddo

  endif

  nx = a(1,1)
  ny = a(2,2)
  nz = a(3,3)

! shifts in the box, modulo 1

  do i = 1,3
    s(i) = lmat(i,1)*mpshift(1) + lmat(i,2)*mpshift(2) + lmat(i,3)*mpshift(3)
  enddo
  if(.not. ldiag) then
    do i = 1,3
      s(i) = s(i) - floor(s(i) + EPS)
      if(abs(s(i)) < EPS) s(i) = ZERO
    enddo
  endif

  sx = s(1)
  sy = s(2)
  sz = s(3)

  return

end subroutine kgrid_box
