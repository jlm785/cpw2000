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

!>  Copies the wave-function coeficients from one k-point to another k-point.
!>  Luttinger-Kohn wavefunction, PR 97, 869 (1954)
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         6 February 2014. 2 November 2025.
!>  \copyright    GNU Public License v2

subroutine psi_convert(neig, nspin,                                      &
     mtxd0, isort0, psi0,   mtxd, isort, psi,                            &
     mxddim, mxdbnd)

! Written 6 February 2014. jlm
! Cleaned and bug mtxd0, 13 April 2019. JLM
! Documentation,February 2020. JLM
! Merged with spin-orbit version. 2 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  mtxd0                           !<  dimension of the hamiltonian (reference)
  complex(REAL64), intent(in)        ::  psi0(nspin*mxddim,nspin*mxdbnd) !<  component j of eigenvector i (reference)
  integer, intent(in)                ::  isort0(mxddim)                  !<  g-vector associated with row/column i of hamiltonian (reference)
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian


! output

  complex(REAL64), intent(out)       ::  psi(nspin*mxddim,nspin*mxdbnd)  !<  component j of eigenvector i (converted)

! local variables

  integer    ::  nmax

! local allocatable arrays

  integer, allocatable               ::  ib0(:)                     !  reverse mapping of isort0

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer ::  i, n


  if(nspin /=1 .and. nspin /= 2) then
    write(6,*)
    write(6,*) '    STOPPED in psi_convert, nspin = ', nspin
    write(6,*)

    STOP

  endif

! The structure of the basis may not be the same, so it is necessary
! to make them compatible.

  nmax = 0
  do i = mtxd,1,-1
    if(isort(i) > nmax) nmax = isort(i)
  enddo
  do i = mtxd0,1,-1
    if(isort0(i) > nmax) nmax = isort0(i)
  enddo

  allocate(ib0(nmax))

  do i=1,nmax
    ib0(i) = 0
  enddo
  do i=1,mtxd0
    ib0(isort0(i)) = i
  enddo

  if(nspin == 1) then

    do n=1,neig
    do i=1,mtxd
      if(ib0(isort(i)) < 1 .or. ib0(isort(i)) > mtxd0) then
        psi(i,n) = C_ZERO
      else
        psi(i,n) = psi0(ib0(isort(i)),n)
      endif
    enddo
    enddo

  else

    do n=1,2*neig
    do i=1,mtxd
      if(ib0(isort(i)) < 1 .or. ib0(isort(i)) > mtxd0) then
        psi(2*i-1,n) = C_ZERO
        psi(2*i  ,n) = C_ZERO
      else
        psi(2*i-1,n) = psi0(2*ib0(isort(i))-1,n)
        psi(2*i  ,n) = psi0(2*ib0(isort(i))  ,n)
      endif
    enddo
    enddo

  endif

  deallocate(ib0)


  return

end subroutine psi_convert
