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

!>  Calculates the matrix <Psi|p|Psi> = <Psi|(1/i) nabla|Psi>.  Complex version
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         January 2013. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine psi_p_psi(mtxd, neig, psi, pmat, rkpt, isort, ng, kgv, nspin, &
    mxddim, mxdbnd, mxdgve)

! Written January 2013. jlm
! Modified (style) 18 February 2014. jlm
! Modified, documentation, 19 February 2020. JLM
! Modified, indentation, lso -> nspin, 4 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands (not counting spin)
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (not counting spin)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions (including spin)
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction
  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  complex(REAL64), intent(in)        ::  psi(nspin*mxddim,nspin*mxdbnd)  !<  wavevector

! output

  complex(REAL64), intent(out)     ::  pmat(3,nspin*mxdbnd,nspin*mxdbnd) !<  <Psi|p|Psi>

! local allocatable arrays

  complex(REAL64), allocatable       ::  p_psi(:,:)                      !  p_j|psi>
  complex(REAL64), allocatable       ::  pmat_j(:,:)                     !  <psi|p_j|psi>

! local variables

  real(REAL64)   ::  qk

! counters

  integer   ::   i, j, k

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  if(nspin /=1 .and. nspin /= 2) then
    write(6,*)
    write(6,*) '    STOPPED in psi_p_psi, nspin = ', nspin
    write(6,*)

    STOP

  endif


  allocate(p_psi(nspin*mtxd,neig))
  allocate(pmat_j(neig,neig))

  do j = 1,3

    do i = 1,neig

      if(nspin == 2) then

        do k = 1,mtxd
          qk = rkpt(j) + UM*(kgv(j,isort(k)))
          p_psi(2*k-1,i) = qk*psi(2*k-1,i)
          p_psi(2*k  ,i) = qk*psi(2*k  ,i)
        enddo

      else

        do k = 1,mtxd
          qk = rkpt(j) + UM*(kgv(j,isort(k)))
          p_psi(k,i) = qk*psi(k,i)
        enddo

      endif

    enddo

    call zgemm('c','n', neig, neig, nspin*mtxd, C_UM, psi, nspin*mxddim, &
               p_psi, nspin*mtxd, C_ZERO, pmat_j, neig)

    do i = 1,neig
    do k = 1,neig
      pmat(j,k,i) = pmat_j(k,i)
    enddo
    enddo

  enddo

  deallocate(p_psi)
  deallocate(pmat_j)

  return

end subroutine psi_p_psi
