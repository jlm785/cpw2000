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
!>  nvec xvec vectors.
!>  If also=.T. performs the corresponding transformation in hxvec.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         January 17 2008, March 5 2025.
!>  \copyright    GNU Public License v2

recursive subroutine grsch_mm_2_c16(also, xvec, hxvec, mtxd, nvec,       &
      irow, xscl,                                                        &
      mxddim)

! It uses a mixture of block classical Gram-Schmidt (using blas3)
! and local modified Gram-Schmidt. All steps are repeated
! according to the "twice is enough" principle.

! Written January 17 2008. JLM
! Modified (f90, complex) February 13 2014.
! Modified, documentation, January 2020. JLM
! Modified, name, indentation, xnorm(nn+n)=0, March 5 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian

  logical, intent(in)                ::  also                            !<  indicates if hxvec is also transformed
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  nvec                            !<  number of vectors
  real(REAL64), intent(in)           ::  xscl(nvec)                      !<  natural scale (norm) of vector i

! input and output

  complex(REAL64), intent(inout)     ::  xvec(mxddim,nvec)               !<  component j of vector i
  complex(REAL64), intent(inout)     ::  hxvec(mxddim,nvec)              !<  component j of companion vector i
  integer, intent(inout)             ::  irow(nvec)                      !<  if irow(i)=0 the vector is linearly dependent on the others and the corresponding xvec is zero


! allocatable arrays

  complex(REAL64), allocatable       ::  tvec(:,:)
  real(REAL64), allocatable          ::  xnorm(:)

! local variables

  real(REAL64)   ::  xsum
  logical        ::  lrepeat

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::   n, nn, nn2

! functions

  complex(REAL64), external  ::  zdotc


  if(nvec <= 10) then

!   calls modified gram-schmidt

    call grsch_mod_c16(also, xvec, hxvec, mtxd, nvec, irow, xscl,        &
                    mxddim)

  else

    nn = nvec/2
    nn2 = nvec - nn

    allocate(tvec(nn,nn2))
    allocate(xnorm(nvec))

!   calls the twin for first half of vectors

    call grsch_mm_c16(also, xvec, hxvec, mtxd, nn, irow, xscl,           &
                      mxddim)

!   orthogonalizes second half to first half

    do n = nn+1,nvec
      xnorm(n) = real( zdotc(mtxd, xvec(1,n), 1, xvec(1,n), 1), REAL64)
    enddo

    call zgemm('c', 'n', nn, nn2, mtxd, C_UM, xvec, mxddim,              &
                xvec(1, nn+1), mxddim, C_ZERO, tvec, nn)

    call zgemm('n', 'n', mtxd, nn2, nn, -C_UM, xvec, mxddim,             &
                tvec, nn, C_UM, xvec(1, nn+1), mxddim)

    if(also) then
      call zgemm('n', 'n', mtxd, nn2, nn, -C_UM, hxvec, mxddim,          &
                  tvec, nn, C_UM, hxvec(1, nn+1), mxddim)
    endif

!   repeat for stability

    lrepeat = .FALSE.
    do n=1,nn2
      xsum = real( zdotc(nn, tvec(1,n), 1, tvec(1,n), 1), REAL64)
      if( xsum > 0.9999*xnorm(nn+n)) lrepeat = .TRUE.
    enddo

    if(lrepeat) then

      call zgemm('c', 'n', nn, nn2, mtxd, C_UM, xvec, mxddim,            &
                  xvec(1, nn+1), mxddim, C_ZERO, tvec, nn)

      call zgemm('n', 'n', mtxd, nn2, nn, -C_UM, xvec, mxddim,           &
                  tvec, nn, C_UM, xvec(1, nn+1), mxddim)

      if(also) then
        call zgemm('n', 'n', mtxd, nn2, nn, -C_UM, hxvec, mxddim,        &
                  tvec, nn, C_UM, hxvec(1, nn+1), mxddim)
      endif

    endif

    deallocate(tvec)
    deallocate(xnorm)

!   calls the twin for second half of vectors

     call grsch_mm_c16(also, xvec(1,nn+1), hxvec(1,nn+1), mtxd, nn2,     &
            irow(nn+1), xscl(nn+1),                                      &
            mxddim)

  endif

  return

end subroutine grsch_mm_2_c16
