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

!>  calculates the product of separable non-local
!>  pseudopotential times neig wavevectors.  complex*16 version
!>  works for both spin and non-spin
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         25 july 2002. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine hk_psi_nl_c16(mtxd, neig, psi, hpsi, anlga, xnlkb, nanl,      &
  ladd, mxddim, mxdbnd, mxdanl)

! written 25 july 2002. jlm
! modified for complex variables and f90  18 June 2012
! Modified, documentation, January 2020. JLM
! modified 7 January 2014, style. JLM
! Modified, indentation, 4 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  integer, intent(in)                ::  nanl                            !<  number of projectors
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  Kleinman-Bylander normalization
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector

  logical, intent(in)                ::  ladd                            !<  true: adds to existing hpsi, false: input hpsi is zeroed

! output

  complex(REAL64), intent(inout)     ::  hpsi(mxddim,mxdbnd)             !<  |hpsi> =  V_NL |psi>

! local variables

  complex(REAL64),allocatable ::  dhd(:,:)
  complex(REAL64) ::  zu

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer   ::   i,n

  if(ladd) then
    zu = C_UM
  else
    zu = C_ZERO
  endif

  allocate(dhd(nanl,neig))

  if(nanl > 0) then

!   dhd = < anl | psi >

    call zgemm('c','n', nanl, neig, mtxd, C_UM, anlga, mxddim,           &
               psi, mxddim, C_ZERO, dhd, nanl)

!   dhd := Diag(xnl) dhd

    do n=1,neig
      do i=1,nanl
        dhd(i,n) = xnlkb(i)*dhd(i,n)
      enddo
    enddo

!   | hpsi > := | hpsi> + | anl > dhd = | hpsi> + | anl > Diag(xnl) < anl | psi >

    call zgemm('n', 'n', mtxd, neig, nanl, C_UM, anlga, mxddim,          &
               dhd, nanl, zu, hpsi, mxddim)

  endif

  deallocate(dhd)

  return

end subroutine hk_psi_nl_c16
