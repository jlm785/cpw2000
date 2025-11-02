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

!>  Calculates the variational k.p (Langreth-Kohn) energies and eigenvectors.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         7 February 2014. 2 November 2025.
!>  \copyright    GNU Public License v2

subroutine kdotp_diag(emax, rkpt, neig, nspin,                           &
    rk0, mtxd0, isort0, psi0, h0, dh0drk, d2h0drk2,                      &
    psi, ei, mtxd, isort, qmod, ekpg,                                    &
    ng, kgv, adot,                                                       &
    mxdgve, mxddim, mxdbnd)


! Written 7 February 2014. JLM
! modified, documentation, February 2020. JLM
! modified, indentation, merged with spin-orbit version, 2 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)

  integer, intent(in)                ::  mtxd0                           !<  dimension of the reference hamiltonian
  integer, intent(in)                ::  isort0(mxddim)                  !<  g-vector associated with row/column i of reference hamiltonian
  complex(REAL64), intent(in)        ::  psi0(nspin*mxddim,nspin*mxdbnd) !<  component j of reference eigenvector

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  rk0(3)                          !<  component in lattice coordinates of the reference k-point
  complex(REAL64), intent(in)        ::  h0(nspin*mxdbnd,nspin*mxdbnd)   !<  <Psi|p|Psi>
  complex(REAL64), intent(in)        ::  dh0drk(nspin*mxdbnd,nspin*mxdbnd,3)      !<  d <Psi|V_NL|Psi> d k
  complex(REAL64), intent(in)        ::  d2h0drk2(nspin*mxdbnd,nspin*mxdbnd,3,3)  !<  d^2 <Psi|V_NL|Psi> d k^2

! output

  complex(REAL64), intent(out)       ::  psi(nspin*mxddim,nspin*mxdbnd)  !<  component j of eigenvector i
  real(REAL64), intent(out)          ::  ei(nspin*mxdbnd)                !<  eigenvalues (Hartree)
  integer, intent(out)               ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(out)               ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(out)          ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(out)          ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

! allocatable arrays

  complex(REAL64), allocatable       ::  psi_lk(:,:)                     !  component j of eigenvector i
  complex(REAL64), allocatable       ::  vec(:,:)                        !  overlap matrix S
  complex(REAL64), allocatable       ::  hred(:,:)                       !  reduced hamiltonian

! other variables

  real(REAL64)      ::  vcell, bdot(3,3)

  integer           ::  info

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::  i, j, m, n


  if(nspin /=1 .and. nspin /= 2) then
    write(6,*)
    write(6,*) '    STOPPED in kdotp_diag, nspin = ', nspin
    write(6,*)

    STOP

  endif

  call adot_to_bdot(adot,vcell,bdot)

  call hamilt_struct(emax, rkpt, mtxd, isort, qmod, ekpg, .FALSE.,       &
      ng, kgv, adot,                                                     &
      mxdgve, mxddim)

  allocate(psi_lk(nspin*mxddim,nspin*mxdbnd))

  call psi_convert(neig, nspin, mtxd0, isort0, psi0, mtxd, isort, psi_lk,  &
        mxddim, mxdbnd)

  allocate(hred(nspin*neig,nspin*neig))
  allocate(vec(nspin*neig,nspin*neig))

  do i=1,nspin*neig
  do j=1,nspin*neig
    hred(j,i) = h0(j,i)
    do n=1,3
      hred(j,i) =  hred(j,i) + dh0drk(j,i,n)*(rkpt(n)-rk0(n))
    enddo
    do n=1,3
    do m=1,3
      hred(j,i) =  hred(j,i) + (rkpt(m)-rk0(m))*d2h0drk2(j,i,m,n)*(rkpt(n)-rk0(n)) / 2
    enddo
    enddo
  enddo
  enddo


  call diag_c16(nspin*neig, hred, ei, vec, nspin*neig, info)


  if(info /= 0) then

    write(6,*) '    Stopped in kdotp_diag'
    write(6,*) '    Diagonalization failed with code: ',info

    stop

  endif


  call zgemm('n','n', nspin*mtxd, nspin*neig, nspin*neig,                &
           C_UM, psi_lk, nspin*mxddim,                                   &
           vec, nspin*neig, C_ZERO, psi, nspin*mxddim)

  deallocate(vec)
  deallocate(hred)

  deallocate(psi_lk)

  return

end subroutine kdotp_diag
