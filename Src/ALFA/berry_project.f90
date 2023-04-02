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

!>  projects several state into our out of a subspace described by orthogonal
!>  states.  For a single state use berry_project_one.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.06
!>  \date         3 February 2023.
!>  \copyright    GNU Public License v2

subroutine berry_project(inout, psi, phi, mtxd, nband, nstate,           &
      mxddim, mxdbnd)

! based on phonon and graham-Schmidt subroutines.

  implicit none

  integer, parameter             :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  character(len=*)                   ::  inout                           !<  'I' or 'O' projects into our out.  If not recognized projects into.

  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  |psi> (in principle eigen-functions)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  nband                           !<  number of bands in space for orthogonalization
  integer, intent(in)                ::  nstate                          !<  number of states to be projected

! input and output

  complex(REAL64), intent(inout)     ::  phi(mxddim,nstate)              !<  wavefunction to be projected.

! local allocatable variables

  complex(REAL64), allocatable       ::  zsum(:,:)                       !  <psi|phi>

! parameters

  real(REAL64), parameter        :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter     :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter     :: C_UM = cmplx(UM,ZERO,REAL64)


  allocate(zsum(nband,nstate))

  call zgemm('C','N', nband, nstate, mtxd, C_UM, psi, mxddim,            &
                     phi, mxddim, C_ZERO, zsum, nband)

  if(inout(1:1) == 'I' .or. inout(1:1) == 'i') then
    call zgemm('N','N', mtxd, nstate, nband, C_UM, psi, mxddim,          &
                     zsum, nband, C_ZERO, phi, mxddim)
  else
    call zgemm('N','N', mtxd, nstate, nband, -C_UM, psi, mxddim,         &
                     zsum, nband, C_UM, phi, mxddim)
  endif

  deallocate(zsum)

  return

end subroutine berry_project
