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
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         February 18, 2014. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine hk_psi_kin_c16(mtxd, neig, psi, hpsi, ekpg, ladd, nspin,      &
    mxddim, mxdbnd)

! Written February 18, 2014, from hk_psi_c16.   jlm
! See that file for historical record.
! Modified, documentation, 26 January 2020. JLM
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! Added spin polarization option, 4 November 2025. JLM


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves (not counting spin)
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands (not counting spin)

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size, not counting spin)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions (including spin)
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  logical, intent(in)                ::  ladd                            !<  true: adds to existing hpsi, false: input hpsi is zeroed

  complex(REAL64), intent(in)        ::  psi(nspin*mxddim,nspin*mxdbnd)  !<  wavevector

! input and output

  complex(REAL64), intent(inout)     ::  hpsi(nspin*mxddim,nspin*mxdbnd) !<  |hpsi> =  K |psi>

! counters

  integer   ::   i, n


  if(nspin /=1 .and. nspin /= 2) then
    write(6,*)
    write(6,*) '    STOPPED in hk_psi_kin_c16, nspin = ', nspin
    write(6,*)

    STOP

  endif

  if(nspin == 1) then                                                    !   no spin


    if(ladd) then

      do n = 1,neig

!$omp parallel do default(shared) private(i)
        do i = 1,mtxd
          hpsi(i,n) = hpsi(i,n) + ekpg(i)*psi(i,n)
        enddo
!$omp end parallel do

      enddo

    else

      do n = 1,neig

!$omp parallel do default(shared) private(i)
        do i = 1,mtxd
          hpsi(i,n) = ekpg(i)*psi(i,n)
        enddo
!$omp end parallel do

      enddo

    endif

  else                                                                      !   with spin

    if(ladd) then

      do n = 1,neig

!$omp parallel do default(shared) private(i)
        do i = 1,mtxd
          hpsi(2*i-1,n) = hpsi(2*i-1,n) + ekpg(i)*psi(2*i-1,n)
          hpsi(2*i  ,n) = hpsi(2*i  ,n) + ekpg(i)*psi(2*i  ,n)
        enddo
!$omp end parallel do

      enddo

    else

      do n = 1,neig

!$omp parallel do default(shared) private(i)
        do i = 1,mtxd
          hpsi(2*i-1,n) = ekpg(i)*psi(2*i-1,n)
          hpsi(2*i  ,n) = ekpg(i)*psi(2*i  ,n)
        enddo
!$omp end parallel do

      enddo

    endif

  endif                                                                      !   with spin

  return

end subroutine hk_psi_kin_c16
