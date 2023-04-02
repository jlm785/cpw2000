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

!>  projects one state into our out of a subspace described by orthogonal
!>  states.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.06
!>  \date         18 January 2023.
!>  \copyright    GNU Public License v2

subroutine berry_project_one(inout, psi, phi, mtxd, nband,               &
      mxddim, mxdbnd)

! based on phonon subroutines

  implicit none

  integer, parameter             :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  character(len=*)                   ::  inout                           !<  'I' or 'O' projects into our out.  If not recognized projects into.

  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  |psi> (in principle eigen-functions)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  nband                           !<  number of bands

! input and output

  complex(REAL64), intent(inout)     ::  phi(mxddim)                     !<  wavefunction to be projected.

! local allocatable variables

  complex(REAL64), allocatable       ::  zsum(:)                         !  <psi|phi>

! parameters

  real(REAL64), parameter        :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter     :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter     :: C_UM = cmplx(UM,ZERO,REAL64)


  allocate(zsum(nband))

  call zgemv('C', mtxd, nband, C_UM, psi, mxddim, phi, 1, C_ZERO, zsum, 1)

  if(inout(1:1) == 'I' .or. inout(1:1) == 'i') then
    call zgemv('N', mtxd, nband, C_UM, psi, mxddim, zsum, 1, C_ZERO, phi, 1)
  else
    call zgemv('N', mtxd, nband,-C_UM, psi, mxddim, zsum, 1, C_UM, phi, 1)
  endif

  deallocate(zsum)

  return

end subroutine berry_project_one
