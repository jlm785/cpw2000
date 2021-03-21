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

!>  Calculates the product of the kinetic energy operator times
!>  neig wavevectors.  complex version

subroutine hk_psi_kin_c16(mtxd, neig, psi, hpsi, ekpg, ladd,             &
    mxddim, mxdbnd)

! Written February 18, 2014, from hk_psi_c16.   jlm
! See that file for historical record.
! Modified, documentation, 26 January 2020. JLM
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! copyright INESC-MN/Jose Luis Martins

! version 4.99

  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
  logical, intent(in)                ::  ladd                            !<  true: adds to existing hpsi, false: input hpsi is zeroed

  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector

! input and output

  complex(REAL64), intent(inout)     ::  hpsi(mxddim,mxdbnd)             !<  |hpsi> =  K |psi>

! counters

  integer   ::   i, n


  if(ladd) then

    do n=1,neig

!$omp parallel do default(shared) private(i)
      do i=1,mtxd
        hpsi(i,n) = hpsi(i,n) + ekpg(i)*psi(i,n)
      enddo
!$omp end parallel do

    enddo

  else

    do n=1,neig

!$omp parallel do default(shared) private(i)
      do i=1,mtxd
        hpsi(i,n) = ekpg(i)*psi(i,n)
      enddo
!$omp end parallel do

    enddo

  endif

  return
end subroutine hk_psi_kin_c16
