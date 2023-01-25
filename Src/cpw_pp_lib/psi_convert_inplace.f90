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

!>  Transforms the coeficients of one k-point to the representation
!>  of another k-point.
!>  Luttinger-Kohn wavefunction, PR 97, 869 (1954).
!>
!>  Trades memory savings for longer execution.
!>  Use psi_convert if memory is not an issue.
!>
!>  Uses insertion sorting because arrays are already mostly ordered.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         13 December 2022.
!>  \copyright    GNU Public License v2

subroutine psi_convert_inplace(neig,mtxd_in,isort_in,mtxd_out,isort_out,psi,   &
     mxddim,mxdbnd)

! Adapted for psi_convert. 13 December 2022. JLM
! Bug in/out squashed, 14 December 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  integer, intent(in)                ::  mtxd_in                         !<  dimension of the hamiltonian (input psi)
  integer, intent(in)                ::  isort_in(mxddim)                !<  g-vector associated with row/column i of hamiltonian (input psi)
  integer, intent(in)                ::  mtxd_out                        !<  dimension of the hamiltonian (output psi)
  integer, intent(in)                ::  isort_out(mxddim)               !<  g-vector associated with row/column i of hamiltonian (output psi)


! input and output

  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i

! local variables

  integer    ::  nmax
  integer    ::  ibtmp
  integer    ::  insert
  integer    ::  nfill

! local allocatable arrays

  integer, allocatable               ::  iback(:)                        !  reverse mapping of isort
  integer, allocatable               ::  irank(:)                        !  target order
  integer, allocatable               ::  ifill(:)                        !  there is information for that component

  complex(REAL64), allocatable       ::  ctmp(:)                         !  temporary storage of

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer ::  i, j, jj


! The structure of the basis may not be the same, so it is necessary
! to make them compatible.

  nmax = 0
  do i = mtxd_out,1,-1
    if(isort_out(i) > nmax) nmax = isort_out(i)
  enddo
  do i = mtxd_in,1,-1
    if(isort_in(i) > nmax) nmax = isort_in(i)
  enddo

  allocate(iback(nmax))

  iback(:) = 0
  do i=1,mtxd_out
    iback(isort_out(i)) = i
  enddo

! fills the rank index

  allocate(irank(mtxd_in))
  do j = 1,mtxd_in
    if(iback(isort_in(j)) /= 0) then
      irank(j) = iback(isort_in(j))
    else
      irank(j) = nmax+j
    endif
  enddo

  iback(:) = 0
  do i=1,mtxd_in
    iback(isort_in(i)) = i
  enddo

  allocate(ifill(mtxd_out))

  ifill(:) = 0
  nfill = 0
  do j = 1,mtxd_out
    if(iback(isort_out(j)) == 0) then
      nfill = nfill + 1
      ifill(j) = nfill
    endif
  enddo

  deallocate(iback)

! first card should be valid.  This should almost never occur.

  if(irank(1) > nmax) then
    do j = 2,mtxd_in
      if(irank(j) > nmax) then
        psi(1,:) = psi(j,:)
        irank(1) = irank(j)
        irank(j) = nmax+j

        exit

      endif
    enddo
  endif

! now insert sort valid cards

  allocate(ctmp(neig))

  jj = 1
  do j = 2,mtxd_in
    if(irank(j) /= 0) then

      jj = jj + 1
      if(irank(j) < irank(jj-1)) then
        ctmp(:) = psi(j,:)
        ibtmp = irank(j)
        insert = 1
        do i = jj-1,1,-1
          insert = i
          if(irank(i) <= ibtmp) then
            insert = insert + 1

            exit

          else
            psi(i+1,:) = psi(i,:)
            irank(i+1) = irank(i)
          endif
        enddo
        psi(insert,:) = ctmp(:)
        irank(insert) = ibtmp
      endif

    endif
  enddo

! they are ordered but not in place

  jj = nfill
  do j = mtxd_out,1,-1
    if(ifill(j) /= 0) then
      psi(j,:) = C_ZERO
      jj = jj - 1

      if(jj == 0) exit

    else
      psi(j,:) = psi(j-jj,:)
    endif
  enddo

  deallocate(ifill)
  deallocate(ctmp)
  deallocate(irank)


  return

end subroutine psi_convert_inplace
