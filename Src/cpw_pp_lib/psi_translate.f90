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

!>  Translates the wave-functions by a reciprocal lattice vector
  
subroutine psi_translate(neig, kgshift,                                  &
  mtxd, isort, psi, mtxd_tr, isort_tr,                                   &
  kgv,                                                                   &
  mxddim, mxdbnd, mxdgve)
  
! Written 17 January 2021. jlm
! copyright  Jose Luis Martins/INESC-MN
  
! version 4.99
  
  implicit none
  
  integer, parameter          :: REAL64 = selected_real_kind(12)
  
  
! input
  
  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  
  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  integer, intent(in)                ::  kgshift(3)                      !<  reciprocal lattice vector shift between k-vectors

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian (reference)
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian (reference)

  integer, intent(in)                ::  mtxd_tr                         !<  dimension of the hamiltonian
  integer, intent(in)                ::  isort_tr(mxddim)                !<  g-vector associated with row/column i of hamiltonian
  
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  
  
! input and output
  
  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (original and converted)
  
! local variables
  
  integer    ::  n1m, n2m, n3m
  integer    ::  kgtmp(3)
  integer    ::  istmp
  
! local allocatable arrays
  
  integer, allocatable               ::  isofkg(:,:,:)                   !  reverse mapping of kgv(isort)
  complex(REAL64), allocatable       ::  psi0(:)                         !  local copy of psi for eigenvector n
  
! parameters
  
  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  
! counters
  
  integer ::  i, n
  
  
! constructs the inverse mapping
  
  n1m = 0
  n2m = 0
  n3m = 0
  do i = 1,mtxd
    if(abs(kgv(1,isort(i))) > n1m) n1m = abs(kgv(1,isort(i)))
    if(abs(kgv(2,isort(i))) > n2m) n2m = abs(kgv(2,isort(i)))
    if(abs(kgv(3,isort(i))) > n3m) n3m = abs(kgv(3,isort(i)))
  enddo
  
  allocate(isofkg(-n1m:n1m,-n2m:n2m,-n3m:n3m))

  isofkg(:,:,:) = 0
  do i = 1,mtxd
    isofkg( kgv(1,isort(i)), kgv(2,isort(i)), kgv(3,isort(i)) ) = i
  enddo
  
  allocate(psi0(mxddim))
  
  do n=1,neig
    psi0(:) = psi(:,n)

    do i=1,mtxd_tr
      kgtmp(1) = kgv(1,isort_tr(i)) + kgshift(1)
      kgtmp(2) = kgv(2,isort_tr(i)) + kgshift(2)
      kgtmp(3) = kgv(3,isort_tr(i)) + kgshift(3)
      if(abs(kgtmp(1)) > n1m .or. abs(kgtmp(2)) > n2m .or. abs(kgtmp(3)) > n3m) then
        psi(i,n) = C_ZERO
      else
        istmp = isofkg(kgtmp(1),kgtmp(2),kgtmp(3))
        if(istmp == 0 ) then
          psi(i,n) = C_ZERO
        else
          psi(i,n) = psi0(istmp)
        endif
      endif
    enddo

  enddo

  deallocate(isofkg)
  deallocate(psi0)

  return
  end subroutine psi_translate
