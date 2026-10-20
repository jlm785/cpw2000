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

!>  Hamiltonian in a non-orthogonal atomic orbital basis whose
!>  eigenvalues are the plane-wave energies of the plane-wave states
!>  that are well represented in the atomic orbital space
!>  (modified tight-binding by projection).
!>
!   In the orthogonalized basis S^(-1/2) the hamiltonian is
!   sum_p e_p |phi_p><phi_p| + Q H_ao Q, where the phi_p are the
!   symmetrically (Lowdin) orthonormalized projections of the selected
!   plane-wave states and Q projects on the rest of the atomic orbital
!   space.  The selected plane-wave energies are reproduced exactly.
!
!   Plane-wave levels (degenerate states) are taken in increasing
!   energy, and a level is selected if the average norm of its states
!   in the atomic orbital space is at least xnmin, and if the
!   projections remain well conditioned (smallest eigenvalue of
!   their overlap at least XLAMIN).  The last level is not used,
!   as the number of plane-wave states may have cut a degenerate level.
!
!   It can be used with spinors (spin-orbit) with nbas = 2*norb.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         20 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_proj_hamilt(ipr, xnmin, nbas, neig, ei, bmat,          &
    hao, s, hmtb, nsel,                                                  &
    mxdbas, mxdbnd)

! Written 20 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbas                          !<  array dimension of the atomic orbital basis
  integer, intent(in)                ::  mxdbnd                          !<  array dimension of plane-wave states

  integer, intent(in)                ::  ipr                             !<  print control
  real(REAL64), intent(in)           ::  xnmin                           !<  minimum average norm in the atomic orbital space of a selected level

  integer, intent(in)                ::  nbas                            !<  number of atomic orbital basis functions
  integer, intent(in)                ::  neig                            !<  number of plane-wave states
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  plane-wave energies, increasing order (hartree)
  complex(REAL64), intent(in)        ::  bmat(mxdbas,mxdbnd)             !<  <ao_i|pw_p>

  complex(REAL64), intent(in)        ::  hao(nbas,nbas)                  !<  hamiltonian in the atomic orbital basis
  complex(REAL64), intent(in)        ::  s(nbas,nbas)                    !<  overlap of the atomic orbital basis

! output

  complex(REAL64), intent(out)       ::  hmtb(nbas,nbas)                 !<  modified tight-binding hamiltonian in the atomic orbital basis
  integer, intent(out)               ::  nsel                            !<  number of selected plane-wave states

! local allocatable arrays

  complex(REAL64), allocatable       ::  s12(:,:), s12_inv(:,:)          !  S^(-1/2), S^(1/2)
  complex(REAL64), allocatable       ::  wrk(:,:)
  real(REAL64), allocatable          ::  ev(:)
  complex(REAL64), allocatable       ::  dmat(:,:)                       !  projections in the orthogonalized basis
  complex(REAL64), allocatable       ::  dsel(:,:)                       !  selected projections
  real(REAL64), allocatable          ::  esel(:)                         !  selected energies
  complex(REAL64), allocatable       ::  gram(:,:), wgram(:,:)           !  overlap of selected projections and eigenvectors
  real(REAL64), allocatable          ::  xlam(:)                         !  eigenvalues of the overlap
  complex(REAL64), allocatable       ::  dorth(:,:)                      !  orthonormalized projections
  complex(REAL64), allocatable       ::  horth(:,:)                      !  hamiltonians in the orthogonalized basis
  complex(REAL64), allocatable       ::  hc(:,:), vc(:,:)                !  hamiltonian in the complement and eigenvectors
  real(REAL64), allocatable          ::  ec(:)                           !  eigenvalues in the complement

! local variables

  integer           ::  nc                                               !  dimension of the complement
  integer           ::  info
  integer           ::  ilev, jlev, ndeg
  real(REAL64)      ::  xnorm
  logical           ::  lok

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  real(REAL64), parameter     ::  TOL = 1.0E-5_REAL64                    !  tolerance for degenerate energies (hartree)
  real(REAL64), parameter     ::  XLAMIN = 0.1_REAL64                    !  smallest eigenvalue of the overlap of selected projections

! counters

  integer   ::  i, j, k


  allocate(s12(nbas,nbas), s12_inv(nbas,nbas), wrk(nbas,nbas), ev(nbas))
  allocate(dmat(nbas,neig), dsel(nbas,nbas), esel(nbas))
  allocate(horth(nbas,nbas))

  call ao_int_GetS12(s, s12, s12_inv, wrk, ev, nbas)

! projections of plane-wave states in the orthogonalized basis

  call zgemm('n', 'n', nbas, neig, nbas, C_UM, s12, nbas,                &
      bmat, mxdbas, C_ZERO, dmat, nbas)

! selects plane-wave levels

  nsel = 0
  ilev = 1
  do while(ilev <= neig)

    jlev = ilev
    do while(jlev < neig)
      if(ei(jlev+1) - ei(ilev) > TOL) exit
      jlev = jlev + 1
    enddo

!   the last level may be incomplete

    if(jlev == neig) exit

    ndeg = jlev - ilev + 1

    xnorm = ZERO
    do j = ilev,jlev
      do i = 1,nbas
        xnorm = xnorm + real(dmat(i,j)*conjg(dmat(i,j)),REAL64)
      enddo
    enddo
    xnorm = xnorm / ndeg

    lok = (xnorm >= xnmin .and. nsel + ndeg <= nbas)

!   the projections must remain linearly independent

    if(lok) then
      do j = 1,ndeg
        dsel(:,nsel+j) = dmat(:,ilev+j-1)
      enddo
      allocate(gram(nsel+ndeg,nsel+ndeg), wgram(nsel+ndeg,nsel+ndeg))
      allocate(xlam(nsel+ndeg))
      call zgemm('c', 'n', nsel+ndeg, nsel+ndeg, nbas, C_UM, dsel, nbas,  &
          dsel, nbas, C_ZERO, gram, nsel+ndeg)
      call diag_c16(nsel+ndeg, gram, xlam, wgram, nsel+ndeg, info)
      if(info /= 0 .or. xlam(1) < XLAMIN) lok = .FALSE.
      deallocate(gram, wgram, xlam)
    endif

    if(lok) then
      do j = 1,ndeg
        esel(nsel+j) = ei(ilev+j-1)
      enddo
      nsel = nsel + ndeg
    endif

    ilev = jlev + 1

  enddo

  if(ipr > 1) then
    write(6,'("   ao_int_proj_hamilt: ",i5," plane-wave states selected ",    &
        &     "out of ",i5)') nsel, neig
  endif

! hamiltonian in the orthogonalized basis

  call zgemm('n', 'n', nbas, nbas, nbas, C_UM, hao, nbas,                &
      s12, nbas, C_ZERO, wrk, nbas)
  call zgemm('n', 'n', nbas, nbas, nbas, C_UM, s12, nbas,                &
      wrk, nbas, C_ZERO, horth, nbas)

  if(nsel > 0) then

!   Lowdin orthonormalization of the selected projections

    allocate(gram(nsel,nsel), wgram(nsel,nsel), xlam(nsel))
    allocate(dorth(nbas,nsel))

    call zgemm('c', 'n', nsel, nsel, nbas, C_UM, dsel, nbas,             &
        dsel, nbas, C_ZERO, gram, nsel)
    call diag_c16(nsel, gram, xlam, wgram, nsel, info)

    if(info /= 0) then
      write(6,*)
      write(6,*) '   STOPPED in ao_int_proj_hamilt: diagonalization of the overlap'
      write(6,*)

      stop

    endif

    do j = 1,nsel
    do i = 1,nsel
      gram(i,j) = C_ZERO
      do k = 1,nsel
        gram(i,j) = gram(i,j) + wgram(i,k)*conjg(wgram(j,k)) / sqrt(xlam(k))
      enddo
    enddo
    enddo

    call zgemm('n', 'n', nbas, nsel, nsel, C_UM, dsel, nbas,             &
        gram, nsel, C_ZERO, dorth, nbas)

!   complement Q = 1 - sum |phi><phi|, its eigenvectors with eigenvalue 1 are the last ones

    nc = nbas - nsel

    call zgemm('n', 'c', nbas, nbas, nsel, -C_UM, dorth, nbas,           &
        dorth, nbas, C_ZERO, wrk, nbas)
    do i = 1,nbas
      wrk(i,i) = wrk(i,i) + C_UM
    enddo
    call diag_c16(nbas, wrk, ev, s12_inv, nbas, info)

!   s12_inv is used as work array, so it is recalculated below

    hmtb(:,:) = C_ZERO

    if(nc > 0) then

      allocate(hc(nc,nc), vc(nc,nc), ec(nc))

      call zgemm('n', 'n', nbas, nc, nbas, C_UM, horth, nbas,            &
          s12_inv(:,nsel+1:nbas), nbas, C_ZERO, wrk, nbas)
      call zgemm('c', 'n', nc, nc, nbas, C_UM, s12_inv(:,nsel+1:nbas),   &
          nbas, wrk, nbas, C_ZERO, hc, nc)
      call diag_c16(nc, hc, ec, vc, nc, info)
      call zgemm('n', 'n', nbas, nc, nc, C_UM, s12_inv(:,nsel+1:nbas),   &
          nbas, vc, nc, C_ZERO, wrk, nbas)

      do k = 1,nc
        do j = 1,nbas
        do i = 1,nbas
          hmtb(i,j) = hmtb(i,j) + ec(k)*wrk(i,k)*conjg(wrk(j,k))
        enddo
        enddo
      enddo

      deallocate(hc, vc, ec)

    endif

    do k = 1,nsel
      do j = 1,nbas
      do i = 1,nbas
        hmtb(i,j) = hmtb(i,j) + esel(k)*dorth(i,k)*conjg(dorth(j,k))
      enddo
      enddo
    enddo

    deallocate(gram, wgram, xlam)
    deallocate(dorth)

!   back to the non-orthogonal basis, H = S^(1/2) H_orth S^(1/2)

    call ao_int_GetS12(s, s12, s12_inv, wrk, ev, nbas)

    call zgemm('n', 'n', nbas, nbas, nbas, C_UM, hmtb, nbas,             &
        s12_inv, nbas, C_ZERO, wrk, nbas)
    call zgemm('n', 'n', nbas, nbas, nbas, C_UM, s12_inv, nbas,          &
        wrk, nbas, C_ZERO, hmtb, nbas)

  else

    hmtb(:,:) = hao(:,:)

  endif

  deallocate(s12, s12_inv, wrk, ev)
  deallocate(dmat, dsel, esel)
  deallocate(horth)

  return

end subroutine ao_int_proj_hamilt
