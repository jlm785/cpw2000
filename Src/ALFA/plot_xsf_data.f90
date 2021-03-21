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

!>  writes function values on a grid in a xsf
!>  input file for later plot of a isosurface vesta or xcrysden

subroutine plot_xsf_data(iotape, lvesta, adot, chd, id, n1,n2,n3)

! Written 1 February 2021. JLM
! Adapted from rho_contour3D

! copyright  J.L.Martins, INESC-MN.

! version 4.99

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  id, n1,n2,n3                    !<  dimensions of array chd

  integer, intent(in)                ::  iotape                          !<  tape number
  logical, intent(in)                ::  lvesta                          !<  indicates that the xsf file follows the vesta orientation

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  complex(REAL64), intent(in)        :: chd(id,n2,n3)                    !<  function in the in the fft real space grid

! local:

  real(REAL64)     ::  avec(3,3), bvec(3,3)
  real(REAL64)     ::  a1(3), a2(3), a3(3)

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  real(REAL64), parameter  :: BOHR = 0.5291772109_REAL64
  
! counters

  integer      ::  i, j, k


! converts stuff

  call adot_to_avec_sym(adot,avec,bvec)

  do k = 1,3
    a1(k) = avec(k,1) * BOHR
    a2(k) = avec(k,2) * BOHR
    a3(k) = avec(k,3) * BOHR
  enddo

! writes preamble

  write(iotape,*)
  write(iotape,*) 'BEGIN_BLOCK_DATAGRID_3D'
  write(iotape,*) '  Charge density'

  write(iotape,*) '    BEGIN_DATAGRID_3D_PSI_RHO'

  write(iotape,*) n1,n2,n3

  write(iotape,'(4x,3f20.10)') ZERO,ZERO,ZERO

  if(lvesta) then

! change order for xsf so that 3 axis is along y...

    write(iotape,'(3(1x,f16.10))') a1(2),a1(3),a1(1)
    write(iotape,'(3(1x,f16.10))') a2(2),a2(3),a2(1)
    write(iotape,'(3(1x,f16.10))') a3(2),a3(3),a3(1)

  else

    write(iotape,'(3(1x,f16.10))') a1(1),a1(2),a1(3)
    write(iotape,'(3(1x,f16.10))') a2(1),a2(2),a2(3)
    write(iotape,'(3(1x,f16.10))') a3(1),a3(2),a3(3)

  endif

  do k = 1,n3
  do j = 1,n2
    write(iotape,'(4x,20f15.10)') (real(chd(i,j,k)),i=1,n1)
  enddo
  enddo

  write(iotape,*) 'END_DATAGRID_3D'
  write(iotape,*) 'END_BLOCK_DATAGRID_3D'

  return
end subroutine plot_xsf_data

