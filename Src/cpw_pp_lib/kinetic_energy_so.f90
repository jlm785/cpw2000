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

!>     Calculates the kinetic energy of the wave-functions
!>     spin-orbit version

       subroutine kinetic_energy_so(neig,mtxd,ekpg,psi,ekpsi,            &
     & mxddim,mxdbnd)


!      written february 18 1990. jlm
!      version 4.0. 15 october 93. jlm
!      modified (f90) 14 January 2014. jlm
!      adapted to spin-orbit, 30 June 2014. JLM
!      Modified, documentation, 6 February 2020. JLM
!      copyright INESC-MN/Jose Luis Martins


!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
       integer, intent(in)                ::  neig                       !<  number of eigenvectors (requested on input, modified by degeneracies on output)
       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       real(REAL64), intent(in)           ::  ekpg(mxddim)               !<  kinetic energy (hartree) of k+g-vector of row/column i

       complex(REAL64), intent(in)        ::  psi(2*mxddim,2*mxdbnd)     !<  component j of eigenvector i

!      output

       real(REAL64), intent(out)          ::  ekpsi(2*mxdbnd)            !<  kinetic energy of eigenvector i. (Hartree)

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64

!      counters

       integer    ::  i, j

         do i=1,2*neig
           ekpsi(i) = ZERO
           do j=1,mtxd
             ekpsi(i) = ekpsi(i) + ekpg(j)*                              &
     &          real(psi(2*j-1,i)*conjg(psi(2*j-1,i)) +                  &
     &               psi(2*j  ,i)*conjg(psi(2*j  ,i)),REAL64)
           enddo
         enddo

       return
       end subroutine kinetic_energy_so
