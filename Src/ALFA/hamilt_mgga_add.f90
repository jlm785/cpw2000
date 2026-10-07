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

!>  Adds to the hamiltonian matrix (from hamilt_kb) and to its
!>  diagonal the meta-GGA generalized Kohn-Sham term
!>
!>    < k+G | -1/2 nabla . vtau nabla | k+G' > = 1/2 (k+G).(k+G') vtau(G-G')
!>
!>  with (k+G).(k+G') = ( |k+G|^2 + |k+G'|^2 - |G-G'|^2 ) / 2.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         6 October 2026.
!>  \copyright    GNU Public License v2

subroutine hamilt_mgga_add(mtxd, mtxds, isort, qmod, vtau,               &
    kgv, phase, conj, inds, kmax, indv, ek,                              &
    hdiag, hamk,                                                         &
    mxdgve, mxdnst, mxdcub, mxddim, mxdsml)

! Written 6 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdsml                          !<  array dimension of the hamiltonian matrix

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian diagonal
  integer, intent(in)                ::  mtxds                           !<  dimension of the hamiltonian matrix
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(in)           ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  complex(REAL64), intent(in)        ::  vtau(mxdnst)                    !<  d (rho eps_xc) / d tau for the prototype g-vector in star j

  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

! input and output

  real(REAL64), intent(inout)        ::  hdiag(mxddim)                   !<  hamiltonian diagonal
  complex(REAL64), intent(inout)     ::  hamk(mxdsml,mxdsml)             !<  hamiltonian matrix of the leading mtxds x mtxds block

! local variables

  real(REAL64)       ::  dotkg                                           !  (k+G_i).(k+G_j)
  complex(REAL64)    ::  hh
  integer   ::  iadd, jadd, kadd
  integer   ::  imx, jmx, kmx
  integer   ::  kxd, kyd, kzd

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer   ::  i, j


! diagonal, 1/2 |k+G|^2 vtau(0)

  do i = 1,mtxd
    hdiag(i) = hdiag(i) + qmod(i)*qmod(i)*real(vtau(1),REAL64) / 2
  enddo

  imx = 2*kmax(1) + 1
  jmx = 2*kmax(2) + 1
  kmx = 2*kmax(3) + 1

  do i = 1,mtxds
    iadd = isort(i)

    hamk(i,i) = hamk(i,i) + qmod(i)*qmod(i)*real(vtau(1),REAL64) / 2

    do j = 1,i-1
      jadd = isort(j)
      kxd = kgv(1,iadd)-kgv(1,jadd)+kmax(1)+1
      if (kxd >= 1 .and. kxd <= imx) then
      kyd = kgv(2,iadd)-kgv(2,jadd)+kmax(2)+1
      if (kyd >= 1 .and. kyd <= jmx) then
      kzd = kgv(3,iadd)-kgv(3,jadd)+kmax(3)+1
      if (kzd >= 1 .and. kzd <= kmx) then

        kadd = ((kxd-1)*jmx+kyd-1)*kmx+kzd
        kadd = indv(kadd)
        if (kadd /= 0) then

          dotkg = (qmod(i)*qmod(i) + qmod(j)*qmod(j) - 2*ek(inds(kadd))) / 2

          hh = vtau(inds(kadd))*conjg(phase(kadd))
          if(conj(kadd) < ZERO) hh = conjg(hh)
          hh = (UM/2)*dotkg*hh

          hamk(i,j) = hamk(i,j) + hh
          hamk(j,i) = hamk(j,i) + conjg(hh)

        endif

      endif
      endif
      endif
    enddo
  enddo

  return

end subroutine hamilt_mgga_add
