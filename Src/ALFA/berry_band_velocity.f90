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

!>  Calculates the band velocity d E_i / d k in a given direction,
!>  or any other Berry quantities in a given direction
!>  with rank 1 (not counting the eigenstate index).
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         21 September 2023.
!>  \copyright    GNU Public License v2

subroutine berry_band_velocity(xk, adot, psidhdkpsi, deidxk,             &
    nlevel, levdeg, leveigs,                                             &
    mxdbnd, mxdlev, mxddeg)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlev                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  real(REAL64), intent(in)           ::  xk(3)                           !<  k-direction in reciprocal lattice coordinates

  integer, intent(in)                ::  nlevel                          !<  number of energy levels
  integer, intent(in)                ::  levdeg(mxdlev)                  !<  degeneragy of level
  integer, intent(in)                ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  complex(REAL64), intent(in)        ::  psidhdkpsi(mxddeg,mxddeg,3,mxdlev)  !<  <psi_n| d H / d k |psi_m>  for example (lattice coordinates)

! output

  real(REAL64), intent(out)          ::  deidxk(mxdbnd)                   !<  d E_i /d xk  in xk direction

! allocatable arrays

  complex(REAL64), allocatable       ::  aoper(:,:)                      !  hermitian operator
  complex(REAL64), allocatable       ::  aeigvec(:,:)                    !  eigenvectors
  real(REAL64), allocatable          ::  aeig(:)                         !  eigenvalues

! local variables

  real(REAL64)      ::  vcell, bdot(3,3)
  real(REAL64)      ::  yk(3), xknorm
  integer           ::  info

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64

! counters

  integer    ::  i, j, n
  integer    ::  nl, nk, mk


! normalizes xk

  call adot_to_bdot(adot, vcell, bdot)

  xknorm = ZERO
  do i = 1,3
  do j = 1,3
    xknorm = xknorm + xk(i)*bdot(i,j)*xk(j)
  enddo
  enddo
  do i = 1,3
    yk(i) = xk(i) / sqrt(xknorm)
  enddo

  allocate(aoper(mxddeg,mxddeg))
  allocate(aeigvec(mxddeg,mxddeg))
  allocate(aeig(mxddeg))

! loop over energy levels

  do nl = 1,nlevel

    if(levdeg(nl) == 1) then

      n = leveigs(nl,1)
      deidxk(n) = ZERO
      do j = 1,3
        deidxk(n) = deidxk(n) + real(psidhdkpsi(1,1,j,nl),REAL64)*yk(j)
      enddo
      write(6,'(i5,f12.4)') n, deidxk(n)

    else

      do nk = 1,levdeg(nl)
      do mk = 1,levdeg(nl)
        aoper(nk,mk) = ZERO
        do j = 1,3
          aoper(nk,mk) = aoper(nk,mk) + psidhdkpsi(nk,mk,j,nl)*yk(j)
        enddo
      enddo
      enddo

      call diag_c16(levdeg(nl), aoper, aeig, aeigvec, mxddeg, info)
      if(info /= 0) then

        write(6,*) "   STOPPED in berry_band_velocity  info = ",info

        stop

      endif

      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        deidxk(n) = aeig(nk)
      enddo

    endif
  enddo


  deallocate(aoper)
  deallocate(aeigvec)
  deallocate(aeig)

  return

end subroutine berry_band_velocity
