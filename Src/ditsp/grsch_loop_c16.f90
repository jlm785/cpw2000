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

!>  Performs a Gram-Schmidt orthogonalization step of
!>  nvec xvec vectors. The first nconv vectors are supposed
!>  to be already orthogonal.
!>  Calls BLAS-3 subroutines for speed.
!>  If also=.T. performs the corresponding transformation in hxvec.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         January 17 2008, March 5 2025.
!>  \copyright    GNU Public License v2


subroutine grsch_loop_c16(also, xvec, hxvec, mtxd, nconv, nvec, irow,    &
                          mxddim)

! It uses a mixture of block classical Gram-Schmidt (using blas3)
! and local modified Gram-Schmidt. All steps are repeated
! according to the "twice is enough" principle.

! Written January 17 2008. jlm
! Modified (f90, complex) February 13 2014.
! Modified, documentation, January 2020. JLM
! Modified, indentation, grsch_mm_c16. March 5 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian

  logical, intent(in)                ::  also                            !<  indicates if hxvec is also transformed
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  nvec                            !<  number of vectors
  integer, intent(in)                ::  nconv                           !<  number of vectors that are already orthogonal

! input and output

  complex(REAL64), intent(inout)     ::  xvec(mxddim,nvec)               !<  component j of vector i
  complex(REAL64), intent(inout)     ::  hxvec(mxddim,nvec)              !<  component j of companion vector i

! output

  integer, intent(out)               ::  irow(nvec)                      !<  if irow(i)=0 the vector is linearly dependent on the others and the corresponding xvec is zero

! allocatable arrays

  real(REAL64), allocatable          ::  xscl(:)
  complex(REAL64), allocatable       ::  tvec(:,:)

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::   n

! functions

  complex(REAL64), external  ::  zdotc


  if( nvec <= nconv ) return


  allocate(xscl(nvec))

  do n = nconv+1,nvec
    irow(n) = 1
  enddo

  do n = nconv+1,nvec
    xscl(n) = real( zdotc(mtxd, xvec(1,n), 1, xvec(1,n), 1), REAL64)
  enddo

! orthogonality to the first nconv vectors

  if(nconv > 0) then

    allocate(tvec(nconv,nvec-nconv))

    call zgemm('c','n', nconv, nvec-nconv, mtxd, C_UM, xvec, mxddim,     &
                xvec(1,nconv+1), mxddim, C_ZERO, tvec, nconv)

    call zgemm('n','n', mtxd, nvec-nconv, nconv,-C_UM, xvec, mxddim,     &
                tvec, nconv, C_UM, xvec(1,nconv+1), mxddim)

    if(also) then
      call zgemm('n','n', mtxd, nvec-nconv, nconv,-C_UM, hxvec, mxddim,  &
                tvec, nconv, C_UM, hxvec(1,nconv+1), mxddim)
    endif

!   repeat for stability

    call zgemm('c','n', nconv, nvec-nconv, mtxd, C_UM, xvec, mxddim,     &
                xvec(1,nconv+1), mxddim, C_ZERO, tvec, nconv)

    call zgemm('n','n', mtxd, nvec-nconv, nconv,-C_UM, xvec, mxddim,     &
                tvec, nconv, C_UM, xvec(1, nconv+1), mxddim)

    if(also) then
      call zgemm('n','n', mtxd, nvec-nconv, nconv,-C_UM, hxvec, mxddim,  &
                tvec, nconv, C_UM, hxvec(1, nconv+1), mxddim)
    endif

    deallocate(tvec)

  endif

! now nconv+1 to nvec

  call grsch_mm_c16(also, xvec(1,nconv+1), hxvec(1,nconv+1), mtxd,       &
      nvec-nconv, irow(nconv+1), xscl(nconv+1),                          &
      mxddim)

  deallocate(xscl)

  return

end subroutine grsch_loop_c16

