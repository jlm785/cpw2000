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

!>  Created the little group associated with the reciprocal space
!>  vector q
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         27 November 2022.
!>  \copyright    GNU Public License v2

subroutine sym_little_group_q(qvec,                                      &
      ntrans_q, mtrx_q, tnp_q, lminus_q, irotm_q, lgam_q,                &
      ntrans, mtrx, tnp)

! Adapted from other sym subroutines. 27 November 2022. JLM


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  qvec(3)                         !<  wave-vector in lattice coordinates

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! output

  integer, intent(out)               ::  ntrans_q                        !<  number of symmetry operations in the little group of q
  integer, intent(out)               ::  mtrx_q(3,3,48)                  !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the little group of q
  real(REAL64), intent(out)          ::  tnp_q(3,48)                     !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the little group of q

  integer, intent(out)               ::  irotm_q                         !<  the rotation sending q -> -q
  logical, intent(out)               ::  lgam_q                          !<  TRUE of q=0
  logical, intent(out)               ::  lminus_q                        !<  if TRUE there is an inversion symmetry

! local variables

  real(REAL64)         ::  qrot(3)
  real(REAL64)         ::  distsq, dif(3)
  integer              ::  isum

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-8_REAL64                        !  tolerance

! counters

  integer    ::  i, j, n


  distsq = ZERO
  do i = 1,3
    distsq = distsq + qvec(i)*qvec(i)
  enddo

  if(distsq < 5*EPS) then

!   gamma point

    lgam_q = .TRUE.
    ntrans_q = ntrans

    do n = 1,ntrans
      do i = 1,3
        do j = 1,3
          mtrx_q(i,j,n) = mtrx(i,j,n)
        enddo
        tnp_q(i,n) = tnp(i,n)
      enddo
    enddo

  else

!   other points

    lgam_q = .FALSE.
    ntrans_q = 0

    do n = 1,ntrans

      do i = 1,3
        qrot(i) = ZERO
        do j = 1,3
          qrot(i) = qrot(i) + mtrx(i,j,n)*qvec(j)
        enddo
      enddo
      distsq = ZERO
      do i = 1,3
        dif(i) = qrot(i) - qvec(i)
        dif(i) = dif(i) - nint(dif(i))
        distsq = distsq + dif(i)*dif(i)
      enddo

      if(distsq < EPS) then
        ntrans_q = ntrans_q + 1
        do i = 1,3
          do j = 1,3
            mtrx_q(i,j,ntrans_q) = mtrx(i,j,n)
          enddo
          tnp_q(i,ntrans_q) = tnp(i,n)
        enddo
      endif

    enddo

  endif

! checks for presence of inverse

  irotm_q = 0
  lminus_q = .FALSE.
  do n = 1,ntrans_q

    isum = 0
    do i = 1,3
    do j = 1,3
      if(i == j) then
        isum = isum + iabs(mtrx_q(i,i,n)+1)
      else
        isum = isum + iabs(mtrx_q(i,j,n))
      endif
    enddo
    enddo
    if(isum == 0) then
      lminus_q = .TRUE.
      irotm_q = n

      exit

    endif
  enddo

  return

end subroutine sym_little_group_q


