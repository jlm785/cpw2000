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

!>  calculates the matrix <Psi|V_NL|Psi> for a separable non-local
!>  pseudopotential V_NL for neig wavevectors.  complex version
!>
!>  \author       Jose Luis Martins
!>  \version      5.02
!>  \date         1980s,  12 September 2021, 28 October 2021.
!>  \copyright    GNU Public License v2

subroutine psi_vnl_psi(mtxd, neig, psi, vnl, anlga, xnlkb, nanl,         &
    mxddim, mxdbnd, mxdanl)

! written June 2012. jlm
! Modified 7 January 2014, style. jlm
! Modified (ladd) 6 February 2014. JLM
! Modified, less memory, documentation. 13 January 2020. JLM
! Modified, nanl > 0,  28 October 2021. JLM
! copyright INESC-MN/Jose Luis Martins

! version 5.02

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  wavefunction dimension
  integer, intent(in)                ::  nanl                            !<  number of projectors
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  Kleinman-Bylander normalization

  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevectors

! output

  complex(REAL64), intent(out)       ::  vnl(mxdbnd,mxdbnd)              !<  <Psi|V_NL|Psi>

! local variables

  complex(REAL64),allocatable      ::  dhd(:,:)
  complex(REAL64),allocatable      ::  xdhd(:,:)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer   ::   i,n


  if(nanl > 0) then

    allocate(dhd(nanl,neig))
    allocate(xdhd(nanl,neig))

!   dhd = < anl | psi >

    call zgemm('c', 'n', nanl, neig, mtxd, C_UM, anlga, mxddim, psi,     &
                   mxddim, C_ZERO, dhd, nanl)

!   xdhd := Diag(xnl) dhd

    do n = 1,neig
      do i = 1,nanl
        xdhd(i,n) = xnlkb(i)*dhd(i,n)
      enddo
    enddo

!   <Psi|V_NL|Psi> = < psi | anl > Diag(xnl)  < anl | psi >

    call zgemm('c', 'n', neig, neig, nanl, C_UM, dhd, nanl, xdhd, nanl,  &
                   C_ZERO, vnl, mxdbnd)

    deallocate(dhd)
    deallocate(xdhd)

  else

    do n = 1,neig
    do i = 1,neig
      vnl(i,n) = C_ZERO
    enddo
    enddo

  endif

  return

end subroutine psi_vnl_psi
