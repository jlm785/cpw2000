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

!>   Given |bas> and H|bas> where |bas> is orthonormal
!>   diagonalizes <bas|H|bas>.
!>
!>  \author       Jose Luis Martins
!>  \version      5.09
!>  \date         April 26 2019. 10 December 2023.
!>  \copyright    GNU Public License v2

subroutine diag_psi_h_psi_c16(mtxd, neig, eg, psi, hpsi, nout, lhpsi,    &
      nbas, bas, hbas, epsdeg,                                              &
      mxddim, mxdbnd, mxdorb)

! Adapted April 26 2019. JLM
! Adapted for 5.0X, 10 December 2023.


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !< array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !< array dimension for number of bands
  integer, intent(in)                ::  mxdorb                          !< array dimension for number of orbitals

  integer, intent(in)                ::  mtxd                            !< wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !< number of wavefunctions

  logical, intent(in)                ::  lhpsi                           !< if true calculates hpsi

  integer, intent(in)                ::  nbas                            !< number of vectors in sub-space basis

  complex(REAL64), intent(in)        ::  bas(mxddim,mxdorb)              !< | bas > sub-space basis set, must be orthonormal on entry
  complex(REAL64), intent(in)        ::  hbas(mxddim,mxdorb)             !< H | bas >

  real(REAL64), intent(in)           ::  epsdeg                          !< two eigenvalues are considered quasi degenerate if |e_i - e_j| < epsdeg

! output

  complex(REAL64), intent(out)       ::  psi(mxddim,mxdbnd)              !< |psi> on input, improved |psi> on output
  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !< H |psi>
  real(REAL64), intent(out)          ::  eg(mxdbnd)                      !< <psi| H |psi>
  integer, intent(out)               ::  nout                            !< number of eigenvectors on output  neig <= nout <= mxdbnd

! local variables

  integer               ::  info

! local allocatable arrays

  complex(REAL64), allocatable       ::  vec(:,:)                        !  eigenvectors of reduced hamiltonian
  complex(REAL64), allocatable       ::  hred(:,:)                       !  reduced hamiltonian
  real(REAL64), allocatable          ::  eig(:)                          !  eigenvalues

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer   ::   n


  allocate(hred(mxdorb,mxdorb))
  allocate(vec(mxdorb,mxdorb))
  allocate(eig(mxdorb))

! psi| H | psi>

  call zgemm('c', 'n', nbas, nbas, mtxd, C_UM, bas, mxddim,              &
       hbas, mxddim, C_ZERO, hred, mxdorb)

! diagonalization

  call diag_c16(nbas, hred, eig, vec, mxdorb, info)

  if(info == 0) then

    nout = neig

    if(mxdbnd > neig) then
      do n = neig + 1,mxdbnd

        if(abs(eig(n) - eig(n-1)) > epsdeg) exit

        nout = n
      enddo
    endif

!   change to plane-wave basis

    call zgemm('n', 'n', mtxd, nout, nbas, C_UM, bas, mxddim,            &
           vec, mxdorb, C_ZERO, psi, mxddim)

    if(lhpsi) then
      call zgemm('n', 'n', mtxd, nout, nbas, C_UM, hbas, mxddim,         &
             vec, mxdorb, C_ZERO, hpsi, mxddim)
    endif

    call dcopy(nout, eig, 1, eg, 1)

  else

    nout =  -1000 + info

  endif

  deallocate(eig)
  deallocate(vec)
  deallocate(hred)

  return

end subroutine diag_psi_h_psi_c16
